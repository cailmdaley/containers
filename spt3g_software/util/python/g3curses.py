"""
Wraps some 3G functionality into the curses library 
"""

import curses, curses.ascii
import time
import traceback
import locale
import types
try:
    from IPython.terminal.interactiveshell import TerminalInteractiveShell
    shell = TerminalInteractiveShell()
    ipython_avail = True
except:
    ipython_avail = False
from datetime import date
from collections import OrderedDict

from spt3g.core import G3Frame, G3FrameType, indexmod

locale.setlocale(locale.LC_ALL,"")
system_year = date.today().year

class EndSession(Exception):
    """
    EndSession Exception cleanly ends a curses session without displaying
    an error message.
    """
    pass

class G3ViewerSubwin(object):
    """
    Represents an individual pane of the viewer, along with the associated
    sub-object. If you live in 1992 and your computer doesn't support unicode,
    you can set current_year=1992 to use ASCII.
    """
    def __init__(self, g3obj, subwindow, width=None, current_year=system_year):
        # Must be associated with a G3ViewerObject
        self.g3obj = g3obj
        self.subwindow = subwindow
        self.startwidth = width
        self.current_year = current_year
        self.scroll_index = 0

    @property
    def index(self):
        return self.g3obj.index

    @index.setter
    def index(self, val):
        self.g3obj.index = val

    @property
    def width(self):
        if self.startwidth is not None:
            return self.startwidth
        return int(self.subwindow.getmaxyx()[1]-self.subwindow.getbegyx()[1])

    def get_object_text(self, index):
        if self.g3obj.can_open(index):
            if self.current_year >= 1993:
                suffix = u'\u27a4'
            else:
                suffix = '>'
        else:
            suffix = u' '
        obj_name = self.g3obj[index]
        obj_text = u'{obj:{leng}.{leng}s}{suffix}\n'.format(
                   obj=obj_name, leng=self.width-2, suffix=suffix)
        return obj_text.encode("utf-8")

class G3ViewerObject:
    """
    Wrapper around an iterable object displayed within a pane. Determines
    whether each constituent item can be further explored.
    """
    def __init__(self, objects):
        self.objects = objects
        self.index = 0

    def __str__(self):
        return str(self.objects)

    def __format__(self, format_mod):
        return str(self).__format__(format_mod)

    def __iter__(self):
        for obj in self.objects:
            yield obj

    def __len__(self):
        return len(self.objects)

    def __getitem__(self, index):
        if self.is_dictlike:
            return G3ViewerObject(self.objects.keys()[index])
        else:
            return G3ViewerObject(self.objects[index])

    @property
    def is_dictlike(self):
        """
        True if the primary viewer object has a values attribute.
        """
        return (hasattr(self.objects, 'values') and
                    (isinstance(self.objects.values, types.BuiltinMethodType) or
                     isinstance(self.objects.values, types.MethodType)))

    @property
    def is_list_of_frames(self):
        """
        True if the primary viewer object is a list of G3Frames.
        """
        try:
            return all([isinstance(entry, G3Frame) for entry in self.objects])
        except:
            return False

    @property
    def can_open_selected(self):
        """
        True if the contents of the selected object can be opened in a new
        window.
        """
        return self.can_open(self.index)

    @property
    def selected_obj_descr(self):
        """
        Returns the description of the selected object.
        """
        if self.is_dictlike:
            obj = self.objects.values()[self.index]
        else:
            obj = self.objects[self.index]
        try:
            return obj.Summary()
        except:
            return str(obj)

    def get_type(self, index):
        """
        Returns the type of the indexed object.
        """
        if self.is_dictlike:
            return type(self.objects.values()[index])
        else:
            return type(self.objects[index])

    @property
    def selected_object(self):
        """
        Returns the object at the selected index as a G3ViewerObject.
        """
        return self.get_object(index=self.index)

    @property
    def selected_object_raw(self):
        """
        Returns the raw object at the selected index.
        """
        return self.get_object_raw(index=self.index)

    def can_open(self, index):
        """
        True if the contents of the indexed object can be expanded in a new
        window.
        """
        if self.is_dictlike:
            return bool(len(self.objects.values()) and
                        hasattr(self.objects.values()[index], '__iter__') and
                        not isinstance(self.objects.values()[index], str))
        else:
            return bool(len(self.objects) and
                        hasattr(self.objects[index], '__iter__') and
                        not isinstance(self.objects[index], str))

    def get_object(self, index):
        """
        Returns the indexed object as a G3ViewerObject.
        """
        if self.is_dictlike:
            return G3ViewerObject(self.objects.values()[index])
        else:
            return G3ViewerObject(self.objects[index])

    def get_object_raw(self, index):
        """
        Returns the raw indexed object.
        """
        if self.is_dictlike:
            return self.objects.values()[index]
        else:
            return self.objects[index]

@indexmod
class ExploreFrame(object):
    """
    Opens each frame in a curses application to allow data exploration.
    """
    def __init__(self, frame_types=[]):
        if not isinstance(frame_types, list):
            log_fatal("frame_types must be a list of G3FrameTypes")
        self.frame_types = frame_types
    def __call__(self, frame):
        if self.frame_types and frame.type not in self.frame_types:
            return
        if frame.type is not G3FrameType.EndProcessing:
            G3Viewer([frame], is_pipeline=True)

class G3Viewer(object):
    """
    Main Curses object explorer. Curses instance will immediately open upon
    instantiation.
    """
    def __init__(self, objects, max_disp_levels=3, current_year=system_year,
                 is_pipeline=False):
        self.max_disp_levels = max_disp_levels
        exit_message = None
        if not isinstance(objects, G3ViewerObject):
            objects = G3ViewerObject(objects)
        self.g3obj = objects
        self.messages = []
        self.subwins = []
        self.obj_stack = []
        self.frame_view = objects.is_list_of_frames
        self._cur_frame = 0
        self.current_year = current_year
        self.is_pipeline = is_pipeline
        # True if the cursor has entered a new subwindow
        self.__redraw_all = True

        try:
            # Init Curses
            self.win = curses.initscr()
            self.win_height = self.win.getmaxyx()[0]-self.win.getbegyx()[1]
            self.win_width = self.win.getmaxyx()[1]-self.win.getbegyx()[1]
            self.msg_win_height = 10
            self.sub_height = self.win_height-self.msg_win_height
            self.sub_width = (self.win_width)//self.max_disp_levels
            self.__init_curses()

            # Create Message Sub-Window:
            self.messagewin = self.win.subpad(10, self.win_width, self.sub_height, 0)

            # Create Types Sub-Window:
            self.typeswin = self.win.subpad(self.sub_height, self.sub_width, 0,
                                        (self.max_disp_levels-1)*self.sub_width)

            # Do the thing
            self.explorer()

        except EndSession:
            pass
        except Exception as e:
            # Store error and print once we leave curses context
            exit_message = traceback.format_exc()
        finally:
            # Clean Up
            self.__deinit_curses()
            if exit_message is not None:
                print(exit_message)

    @property
    def current_frame(self):
        if not self.frame_view:
            # TODO: Raise warning
            return
        return self.g3obj[self._cur_frame].objects

    def __init_curses(self):
        curses.noecho()
        curses.cbreak()
        curses.curs_set(0)
        self.win.keypad(1)

    def __deinit_curses(self):
        self.win.clear()
        curses.nocbreak()
        self.win.keypad(0)
        curses.echo()
        curses.endwin()

    @property
    def current_frame_index(self):
        return self._cur_frame

    @current_frame_index.setter
    def current_frame_index(self, frame_num):
        self._cur_frame = frame_num
        self.subwins[0].g3obj = self.g3obj[frame_num]

    def explorer(self):
        # If the first object is a list of frames, put the frame list in the
        # message window instead of the first object window.
        if self.frame_view:
            iterable = self.g3obj[self.current_frame_index]
        else:
            iterable = self.g3obj
        self.subwins.append(self.__new_win(iterable))
        self.cur_win = 0

        self.__update_buffer()
        while True:
            self.__handle_input()
            self.__update_buffer()

    def __new_types_win(self, obj):
        """
        Right-most window reserved for displaying object properties.
        """
        start = (self.max_disp_levels-1)*self.sub_width
        subw = self.win.subpad(self.sub_height, self.sub_width, 0, start)
        return G3ViewerSubwin(obj, subw, self.sub_width, self.current_year)

    def __new_win(self, obj):
        """
        Open a new subwindow whose content is the expansion of the current
        selected object.
        """
        start = len(self.subwins)*self.sub_width
        subw = self.win.subpad(self.sub_height, self.sub_width, 0, start)
        return G3ViewerSubwin(obj, subw, self.sub_width, self.current_year)

    def __delete_trailing_subwindows(self, index):
        """
        Destroy all subwindows after `index`.
        """
        while index+1 < len(self.subwins):
            old_subw = self.subwins.pop()
            old_subw.subwindow.clear()
            old_subw.subwindow.refresh()

    def __handle_input(self):
        c = self.win.getch()
        self.__redraw_all = False
        subw = self.subwins[self.cur_win]
        if c in [ord('j'), ord('s'), curses.KEY_DOWN]:
            # Go up an item.
            if (subw.index != len(subw.g3obj)-1):
                subw.index += 1
                if subw.index-subw.scroll_index > self.sub_height-2:
                    subw.scroll_index += 1
                    self.__redraw_all = True
                else:
                    # Only change emphasis of current and previous types.
                    type_str = "| {0:.{win_width}s}".format(
                                subw.g3obj.get_type(subw.index),
                                win_width=self.sub_width-2)
                    self.typeswin.addstr(subw.index-subw.scroll_index, 0,
                                         type_str, curses.A_STANDOUT)
                    type_str = "| {0:.{win_width}s}".format(
                                subw.g3obj.get_type(subw.index-1),
                                win_width=self.sub_width-2)
                    self.typeswin.addstr(subw.index-subw.scroll_index-1, 0,
                                         type_str)
                self.__delete_trailing_subwindows(self.cur_win)
        elif c in [ord('g')]:
            # Go to top
            subw.index = 0
            subw.scroll_index = 0
            self.__redraw_all = True
        elif curses.ascii.unctrl(c) in ['^U']:
            # Go up one page
            subw.index = max(subw.index-self.sub_height + 2, 0)
            subw.scroll_index = subw.index
            self.__redraw_all = True
        elif c in [ord('k'), ord('w'), curses.KEY_UP]:
            # Go down an item.
            if (subw.index != 0):
                subw.index -= 1
                if subw.scroll_index == subw.index+1:
                    subw.scroll_index -= 1
                    self.__redraw_all = True
                else:
                    self.__delete_trailing_subwindows(self.cur_win)
                    type_str = "| {0:.{win_width}s}".format(
                                subw.g3obj.get_type(subw.index),
                                win_width=self.sub_width-2)
                    self.typeswin.addstr(subw.index-subw.scroll_index, 0,
                                         type_str, curses.A_STANDOUT)
                    type_str = "| {0:.{win_width}s}".format(
                                subw.g3obj.get_type(subw.index+1),
                                win_width=self.sub_width-2)
                    self.typeswin.addstr(subw.index-subw.scroll_index+1, 0,
                                         type_str)
        elif c in [ord('G')]:
            # Go to bottom
            subw.index = len(subw.g3obj) - 1
            subw.scroll_index = max(subw.index - self.sub_height + 2, 0)
            self.__redraw_all = True
        elif curses.ascii.unctrl(c) in ['^D']:
            # Go down one page
            subw.index = min(subw.index+self.sub_height-2,
                             len(subw.g3obj)-1)
            subw.scroll_index = max(subw.index - self.sub_height + 2, 0)
            self.__redraw_all = True
        elif c in [ord('q')]:
            raise EndSession()
        elif c in [curses.KEY_ENTER, ord('l'), ord('d'), curses.KEY_RIGHT, 10]:
            # If selected item is expandable, open in new subwindow.
            if subw.g3obj.can_open_selected:
                self.__delete_trailing_subwindows(self.cur_win)
                self.__open_in_subwin(subw.g3obj.selected_object)
                self.__redraw_all = True
        elif c in [ord('h'), ord('a'), curses.KEY_LEFT]:
            # If we are not at the top level, go back to parent object.
            if self.cur_win > 0:
                self.cur_win -= 1
            else:
                # Since we will do a full redraw, all we need to do is change
                # the object associated with each subwindow.
                self.__pop_obj_stack()
            self.__redraw_all = True
        elif c in [ord('L'), ord('D'), curses.KEY_SRIGHT]:
            # If we are looking at a list of G3Frames, open the next frame
            if (self.frame_view and
                    self.current_frame_index < len(self.g3obj)-1):
                self.current_frame_index += 1
                self.cur_win = 0
                self.__delete_trailing_subwindows(0)
                self.__redraw_all = True
        elif c in [ord('H'), ord('A'), curses.KEY_SLEFT]:
            # If we are looking at a list of G3Frames, open the previous frame
            if (self.frame_view and self.current_frame_index > 0):
                self.current_frame_index -= 1
                self.cur_win = 0
                self.__delete_trailing_subwindows(0)
                self.__redraw_all = True
        elif c in [ord('o')]:
            # If iPython is available, open the selected object in an
            # interpreter
            if ipython_avail:
                self.__deinit_curses()
                shell.push({
                    'selected': subw.g3obj.selected_object_raw,
                    'current_frame': self.current_frame
                    })
                shell.mainloop(display_banner="""
The current frame (if applicable) is stored in `current_frame`.
The selected object is stored in `selected`.
""")
                self.__init_curses()
                self.__redraw_all = True
        elif c in [ord('?')]:
            self.help_menu()
            self.__redraw_all = True

    def __open_in_subwin(self, objects):
        """
        Opens selected object in right-most subwindow.
        If cuw_win+1==max_disp_levels, will shift all windows left
        and display the new object in the right-most frame. Right-most window
        is reserved for key properties.
        """
        if len(self.subwins)<self.max_disp_levels-1:
            self.subwins.append(self.__new_win(objects))
            self.cur_win += 1
        else:
            self.obj_stack.append(self.subwins[0].g3obj)
            for i in range(self.max_disp_levels-2):
                self.subwins[i].g3obj = self.subwins[i+1].g3obj
            self.subwins[self.max_disp_levels-2].g3obj = objects
        self.cur_win = len(self.subwins)-1

    def __pop_obj_stack(self):
        """
        Call if trying to move left of first frame. If objects exist on stack,
        will pop first to left-most frame and shift rest right, discarding the
        right-most object from the viewer.
        """
        if len(self.obj_stack)>0:
            left_obj = self.obj_stack.pop()
            for i in reversed(range(len(self.subwins)-1)):
                self.subwins[i+1].g3obj = self.subwins[i].g3obj
            self.subwins[0].g3obj = left_obj

    def __update_buffer(self):
        """
        Update the buffer of some or all curses subwindows.
        Setting __redraw_all will cause all windows to be redrawn.
        """
        # Update
        for i, subw in enumerate(self.subwins):
            if self.__redraw_all:
                subw.subwindow.clear()
            elif i<self.cur_win:
                # Windows to our left should never need to be redrawn
                continue

            # Redraw full subwindow
            scroll_index = subw.scroll_index
            iter_len = min(len(subw.g3obj), self.sub_height)
            for j in range(iter_len):
                obj_index = scroll_index + j
                line = subw.get_object_text(obj_index)
                if obj_index == subw.index:
                    if i == self.cur_win:
                        fmt = curses.A_STANDOUT
                    else:
                        fmt = curses.A_BOLD
                else:
                    fmt = curses.A_NORMAL
                subw.subwindow.addstr(j, 0, line, fmt)
                if j+2 == self.sub_height:
                    break
            subw.subwindow.refresh()

        if self.__redraw_all:
            # Update types window for cur_win
            self.typeswin.clear()
            g3obj = self.subwins[self.cur_win].g3obj
            scroll_index = self.subwins[self.cur_win].scroll_index
            iter_len = min(len(g3obj), self.sub_height)
            for i in range(iter_len):
                obj_index = scroll_index + i
                type_str = "| {0:.{win_width}s}".format(
                        g3obj.get_type(obj_index), win_width=self.sub_width-2)
                if obj_index == self.subwins[self.cur_win].index:
                    fmt = curses.A_STANDOUT
                else:
                    fmt = curses.A_NORMAL
                self.typeswin.addstr(i, 0, type_str, fmt)
                if i+2 == self.sub_height:
                    break

            # Write info to message window
            self.messagewin.clear()
            if self.frame_view:
                self.__write_frame_list()
                self.__write_frame_props()
            else:
                message = self.subwins[self.cur_win].g3obj.selected_obj_descr
                self.messagewin.addstr(message.encode("utf-8"))

            self.messagewin.addstr(self.msg_win_height-1, self.win_width-20,
                                   "Press '?' for Help")
            self.msg_win_height = 10
        self.typeswin.refresh()
        self.messagewin.refresh()

    def __write_frame_list(self):
        """
        Write the G3FrameType letter keys for every frame in the viewer, with
        the current frame highlighted.
        """
        # TODO: Scrolling for case len(g3obj) > message window width.
        if self.is_pipeline:
            self.messagewin.addstr('Current Frame Type: ')
            self.messagewin.addstr(self.current_frame.type.name)
        else:
            self.messagewin.addstr('Frame List: ')
            for i, frame in enumerate(self.g3obj):
                if i == self.current_frame_index:
                    fmt = curses.A_STANDOUT
                else:
                    fmt = curses.A_NORMAL
                self.messagewin.addstr(frame.type.key, fmt)
        self.messagewin.addstr('\n')

    def __write_frame_props(self):
        """
        Writes out the "typical" keys and summaries for the selected frame.
        """
        row = 1
        col = 0
        for key, val in get_typical_frame_info(self.current_frame).items():
            self.messagewin.addstr(row, col, "{0:15.15s}: ".format(key), curses.A_BOLD)
            self.messagewin.addstr(row, col+18, "{0:15.15s}".format(val))
            if row+1 == self.msg_win_height:
                row = 1
                col += 35
            else:
                row += 1

    def help_menu(self):
        self.win.clear()
        self.win.addstr("?                Help (Hit Any Key To Leave)\n")
        self.win.addstr("a/h/LEFT         Previous Window\n")
        self.win.addstr("d/l/RIGHT/ENTER  Expand In New Window\n")
        self.win.addstr("w/k/UP           Navigate Up\n")
        self.win.addstr("s/j/DOWN         Navigate Down\n")
        self.win.addstr("A/H/SHIFT+LEFT   Previous Frame\n")
        self.win.addstr("D/L/SHIFT+RIGHT  Next Frame\n")
        self.win.addstr("o                Open Selected Object In iPython\n")
        if self.is_pipeline:
            self.win.addstr("q                Quit (moves on to next frame)\n")
        else:
            self.win.addstr("q                Quit\n")
        c = self.win.getch()
        self.win.clear()

def get_typical_frame_info(frame):
    """
    Save the summary of the "typical" keys (as defined by the all-knowing
    documentation) to a dictionary for display in the message window.
    """
    typical_keys = {
            G3FrameType.Scan:
                ['StartTime', 'StopTime', 'SourceName', 'BoresightAz',
                 'BoresightEl', 'RawTimestreams_I', 'RawTimesreams_Q',
                 'CalTimestreams', 'TimestreamWeights', 'ScanNumber', 'Flags'],
            G3FrameType.Timepoint:
                ['EventHeader', 'DfMux'],
            G3FrameType.Housekeeping:
                ['DfMuxHousekeeping'],
            G3FrameType.Calibration:
                ['BolometerProperties', 'NominalBolometerProperties',
                 'TimeConst'],
            G3FrameType.Observation:
                ['SourceName', 'ObservationNumber'],
            G3FrameType.Wiring:
                ['WiringMap'],
            G3FrameType.GcpSlow:
                ['array', 'antenna0']
            }

    keys = OrderedDict()

    if frame.type in typical_keys:
        for key in typical_keys[frame.type]:
            try:
                keys[key] = frame[key].Summary()
            except:
                pass

    return keys
