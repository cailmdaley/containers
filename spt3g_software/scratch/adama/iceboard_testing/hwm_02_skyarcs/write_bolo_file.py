arciest = [
('ARG1A.1A.38.Y', 283608, 'spaced stripes') ,
('ARG1A.1A.2.X', 287586, 'nice stripes in two') ,
('ARG2B.2B.24.Y', 289692, 'thin random stripes or blotchy') ,
('ARG1A.1A.22.X', 314730, 'corner wide and thin') ,
('ARG1A.1A.24.Y', 318240, 'bright in two') ,
('ARG1B.1B.39.Y', 321984, 'thin bullseye in one') ,
('ARG1A.1A.5.X', 347139, 'nice wide stripes in one') ,
('ARG2B.2B.38.X', 370071, 'arcs but subtle') ,
('ARG2A.2A.17.Y', 376974, 'meh or nice thin stripes') ,
('ARG1A.1A.19.Y', 412191, 'big bright') ,
('ARG2B.2B.22.Y', 414297, 'too horizontal low contrast') ,
('ARG1B.1B.20.X', 424125, 'bright stripes in two') ,
('ARG1B.1B.26.X', 430092, 'bright big') ,
('ARG2B.2B.2.X', 432315, 'thin arc only in two') ,
#('ARG1A.1A.34.X', 444015, 'medium low contrast') ,
('ARG2A.2A.1.Y', 462267, 'meh and archy') ,
('ARG1B.1B.8.X', 465075, 'stripe in one') ,
('ARG1B.1B.23.X', 469404, 'big and little') ,
('ARG1B.1B.32.X', 492336, 'blotchy or stripes') ,
('ARG1A.1A.36.Y', 495612, 'bright in two') ,
('ARG2B.2B.15.X', 500409, 'thin in one place') ,
('ARG1B.1B.9.X', 501696, 'thin stripes in one light arcs in another') ,
('ARG1A.1A.20.Y', 506610, 'big bright or high freq') ,
('ARG1A.1A.18.Y', 511524, 'bright corner zigazgy in one') ,
('ARG2A.2A.18.Y', 516204, 'thin wisps') ,
('ARG2A.2A.20.Y', 517023, 'blotchy') ,
('ARG1A.1A.4.Y', 519012, 'bright arches') ,
('ARG2B.2B.40.Y', 526149, 'two freq stripes crazy arcs') ,
('ARG1B.1B.29.Y', 547911, 'big not strong one stripe') ,
('ARG2A.2A.26.Y', 571077, 'dual direction high freq') ,
('ARG1A.1A.26.Y', 573417, 'thin wisps thick in one') ,
('ARG1A.1A.20.X', 580788, 'one thick stripe') ,
('ARG1A.1A.13.Y', 605358, 'thin arches great pattern') ,
('ARG2A.2A.13.Y', 615069, 'meh or blotchy') ,
('ARG2A.2A.38.X', 620334, 'pattern above noise in one') ,
('ARG1B.1B.3.X', 633438, 'meh or thin stripes') ,
('ARG2B.2B.17.X', 645723, 'one big stripe') ,
('ARG1A.1A.37.X', 652509, 'blotchy arcy') ,
('ARG2B.2B.19.Y', 672867, 'crossing thin stripes in one') ,
('ARG2B.2B.39.X', 673335, 'great in one') ,
('ARG1A.1A.32.X', 673569, 'stripes in noisy obs') ,
('ARG2A.2A.15.X', 680121, 'meh in one') ,
#('ARG2A.2A.9.Y', 682110, 'meh or blotchy') ,
('ARG2A.2A.27.Y', 686907, 'one bright streak') ,
('ARG1A.1A.16.X', 713817, 'not strong but clear in two') ,
('ARG1A.1A.43.X', 725166, 'squiggly corner high freq') ,
('ARG1B.1B.33.X', 728793, 'too horizontal low contrast') ,
('ARG2B.2B.22.X', 740376, 'big stripes blotchy') ,
('ARG2A.2A.11.X', 768456, 'blotchy') ,
('ARG2A.2A.33.Y', 772785, 'meh in one nice thin and thick') ,
('ARG1B.1B.26.Y', 780858, 'variable spacing') ,
('ARG2A.2A.30.X', 783081, 'high contrast changes pattern') ,
('ARG2B.2B.6.Y', 784134, 'big blobs') ,
('ARG1B.1B.2.Y', 785070, 'bright crazy') ,
('ARG1A.1A.1.X', 804609, 'big blobs in one') ,
('ARG2A.2A.39.X', 838422, 'nice stripes only in one thin and thick') ,
('ARG2B.2B.2.Y', 841815, 'low contrast') ,
('ARG2B.2B.15.Y', 844272, 'bright large') ,
('ARG1B.1B.15.Y', 852579, 'few stripes in one') ,
('ARG2B.2B.14.Y', 876213, 'bright stripe in one') ,
('ARG1B.1B.14.Y', 882297, 'bright corners in two') ,
('ARG2B.2B.1.X', 885573, 'one bright line in one') ,
('ARG1B.1B.40.X', 886743, 'nothing in one big bright in others') ,
('ARG2B.2B.7.Y', 915759, 'occasional streak nice thin') ,
('ARG2A.2A.25.Y', 933543, 'wide in two') ,
#('ARG1A.1A.39.Y', 937755, 'unclear in one') ,
('ARG2B.2B.20.Y', 999414, 'good pattern in one high freq') ,
('ARG1B.1B.16.X', 1133379, 'weak circle in two') ,
]

if False:
    print "bolometer\tchannel"
    for jsquid in range(4):
        nchannel = 1
        for bolo in arciest:
            print bolo[0] + '_' + str(jsquid+1) + '\t0112/1/' + str(jsquid+1) + '/' + str(nchannel)
            nchannel = nchannel + 1
            
if True:
    for jsquid in range(4):
        nchannel = 1
        for bolo in arciest:
            print "%s_%d: {\n    frequency: %.1f\n}\n" % (bolo[0],jsquid+1,bolo[1])
