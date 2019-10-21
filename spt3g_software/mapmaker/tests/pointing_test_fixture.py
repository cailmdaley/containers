from __future__ import print_function
import unittest
import numpy as np
import pickle
from spt3g import core, mapmaker, coordinateutils

un = core.G3Units
projs = coordinateutils.MapProjection

def sloppy_eq(a,b, slop = 1e-1):
    return (np.isnan(a) and np.isnan(b) ) or abs(a-b) < slop

def test_inversion(flat_map):
    for pix in range(np.size(flat_map)):
        ang = flat_map.pixel_to_angle(pix)
        pix_o = flat_map.angle_to_pixel(ang[0], ang[1])
        if (pix != pix_o):
            print('returning false')
            print(pix, pix_o)
            return False
    return True


class TestPointingMethods(unittest.TestCase):
    #@unittest.skip("demonstrating skipping")
    def test_pointing(self):
        pointing_test_cases = pickle.load(open('test_pointing_information.pkl'))
        for d in pointing_test_cases:
            flipped_pix_x = d['n_x'] - int(d['pix_x']) -1
            flipped_pix_y = d['n_y'] - int(d['pix_y']) -1
            if flipped_pix_x < 0 or flipped_pix_y < 0:
                continue

            #import pdb; pdb.set_trace()
            fsm = coordinateutils.FlatSkyMap(d['n_x'], d['n_y'],
                                             d['reso_arcmin']*un.arcmin,
                                             proj = coordinateutils.MapProjection(d['proj']), 
                                             alpha_center = d['ra0']*un.deg, 
                                             delta_center = d['dec0']*un.deg
            ) 
            sw_pix = coordinateutils.pixel_1d_to_pixel_2d(
                fsm.angle_to_pixel(d['ra_in']*un.deg, d['dec_in']*un.deg),
                d['n_x'], d['n_y'])

            
            sw_ang = fsm.pixel_to_angle(  flipped_pix_x , flipped_pix_y)
            sw_ang[0]/= un.deg
            sw_ang[1]/= un.deg

            print(sw_ang)

            #print(flipped_pix_x, flipped_pix_y, sw_pix, d['in_ss_pix'], (sw_pix[0] != -1))
            #print( (bool(d['in_ss_pix']) == (sw_pix[0] != -1) ) )
            self.assertTrue( (bool(d['in_ss_pix']) == (sw_pix[0] != -1) )  or
                             (sloppy_eq(sw_pix[0], 0, 2) and not bool(d['in_ss_pix']) )or
                             (sloppy_eq(sw_pix[1], 0, 2) and not bool(d['in_ss_pix']))
                             
                         )
            if ( d['in_ss_pix']):
                print("Ra", d['ra_out'], sw_ang[0])
                print("dec", d['dec_out'], sw_ang[1])
                print("proj", d['proj'])
                
                
                self.assertTrue( sloppy_eq(flipped_pix_x, sw_pix[1],  2) )
                self.assertTrue( sloppy_eq(flipped_pix_y, sw_pix[0],  2) )
                self.assertTrue( sloppy_eq( d['ra_out'], sw_ang[0] % 360.) )
                self.assertTrue( sloppy_eq( d['dec_out'], sw_ang[1]) )

    #@unittest.skip("demonstrating skipping")
    def test_invert(self):
        delta_center = -55.0* core.G3Units.deg
        alpha_center = 125.5 * core.G3Units.deg
        res = 1.0 * core.G3Units.arcmin
        n_x = 780
        n_y = 780        
        proj = coordinateutils.MapProjection.Proj0
        for proj in [coordinateutils.MapProjection.Proj0,
                     coordinateutils.MapProjection.Proj1,
                     coordinateutils.MapProjection.Proj2,
                     coordinateutils.MapProjection.Proj4,
                     coordinateutils.MapProjection.Proj5 ]:
            out_map = coordinateutils.FlatSkyMap(x_len = n_x, y_len = n_y,
                                          res = res,
                                          delta_center = delta_center,
                                          alpha_center = alpha_center,
                                          proj= proj)        
            self.assertTrue(test_inversion(out_map))
    
            
if __name__ == '__main__':
    unittest.main()
