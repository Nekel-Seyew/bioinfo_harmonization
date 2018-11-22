import math
import mpiAlgo
import random
from statistics import mean, stdev

from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
worldsize = comm.Get_size()

#if rank > 0:
#    while True:
#        data = comm.recv(source=0,tag=1)
#        child = mpiAlgo.multiproc_vertex(data)
#        req = comm.isend(child,dest=0,tag=1)
#        req.wait()


windowSize = 17 #input("Please input window size (odd numbers only)")
freqDict = dict() #dictionary mapping codons to their frequencies
mapDict = dict() #dictionary mapping codons to amino acid
aaFreqDict = dict() #dictionary mapping each amino acid to a list of the frequencies of possible codons for that amino acid
aaMapDict = dict() #dictionary from amino acid to list of codons with frequencies for it (for RRTs)
mapDict = {'TCA': 'S', 'AAT': 'N', 'TGG': 'W', 'GAT': 'D', 'GAA': 'E', 'TTC': 'F', 'CCG': 'P',
           'ACT': 'T', 'GGG': 'G', 'ACG': 'T', 'AGA': 'R', 'TTG': 'L', 'GTC': 'V', 'GCA': 'A',
           'TGA': '*', 'CGT': 'R', 'CAC': 'H', 'CTC': 'L', 'CGA': 'R', 'GCT': 'A', 'ATC': 'I',
           'ATA': 'I', 'TTT': 'F', 'TAA': '*', 'GTG': 'V', 'GCC': 'A', 'GAG': 'E', 'CAT': 'H',
           'AAG': 'K', 'AAA': 'K', 'GCG': 'A', 'TCC': 'S', 'GGC': 'G', 'TCT': 'S', 'CCT': 'P',
           'GTA': 'V', 'AGG': 'R', 'CCA': 'P', 'TAT': 'Y', 'ACC': 'T', 'TCG': 'S', 'ATG': 'M',
           'TTA': 'L', 'TGC': 'C', 'GTT': 'V', 'CTT': 'L', 'CAG': 'Q', 'CCC': 'P', 'ATT': 'I',
           'ACA': 'T', 'AAC': 'N', 'GGT': 'G', 'AGC': 'S', 'CGG': 'R', 'TAG': '*', 'CGC': 'R',
           'AGT': 'S', 'CTA': 'L', 'CAA': 'Q', 'CTG': 'L', 'GGA': 'G', 'TGT': 'C', 'TAC': 'Y',
           'GAC': 'D'}

speciesDict = {
        'Escherichia_coli' : ['Table contains 22512846 CDSs (6729157102 codons)\n',
                              '\n',
                              'TTT\xa022.38\xa0(150622203)\xa0\xa0TCT\xa0\xa08.61\xa0(57928149)\xa0\xa0\xa0TAT\xa016.36\xa0(110090438)\xa0\xa0TGT\xa0\xa05.19\xa0(34921543)\xa0\n',
                              'TTC\xa016.21\xa0(109055063)\xa0\xa0TCC\xa0\xa08.81\xa0(59303334)\xa0\xa0\xa0TAC\xa012.15\xa0(81734382)\xa0\xa0\xa0TGC\xa0\xa06.34\xa0(42659065)\xa0\n',
                              'TTA\xa013.83\xa0(93067307)\xa0\xa0\xa0TCA\xa0\xa07.57\xa0(50962811)\xa0\xa0\xa0TAA\xa0\xa02.03\xa0(13685782)\xa0\xa0\xa0TGA\xa0\xa01.04\xa0(6992446)\xa0\xa0\n',
                              'TTG\xa013.37\xa0(90002242)\xa0\xa0\xa0TCG\xa0\xa08.79\xa0(59125602)\xa0\xa0\xa0TAG\xa0\xa00.25\xa0(1691565)\xa0\xa0\xa0\xa0TGG\xa015.21\xa0(102383498)\n', 
                              '\n',
                              'CTT\xa011.44\xa0(77012997)\xa0\xa0\xa0CCT\xa0\xa07.22\xa0(48555280)\xa0\xa0\xa0CAT\xa012.84\xa0(86374193)\xa0\xa0\xa0CGT\xa020.70\xa0(139269987)\n',
                              'CTC\xa010.92\xa0(73471955)\xa0\xa0\xa0CCC\xa0\xa05.56\xa0(37398058)\xa0\xa0\xa0CAC\xa0\xa09.44\xa0(63551268)\xa0\xa0\xa0CGC\xa021.48\xa0(144508930)\n',
                              'CTA\xa0\xa03.93\xa0(26430058)\xa0\xa0\xa0CCA\xa0\xa08.44\xa0(56798739)\xa0\xa0\xa0CAA\xa015.10\xa0(101596569)\xa0\xa0CGA\xa0\xa03.67\xa0(24720122)\xa0\n',
                              'CTG\xa052.10\xa0(350589421)\xa0\xa0CCG\xa022.65\xa0(152409497)\xa0\xa0CAG\xa029.21\xa0(196583131)\xa0\xa0CGG\xa0\xa05.72\xa0(38486931)\xa0\n',
                              '\n',
                              'ATT\xa030.21\xa0(203300538)\xa0\xa0ACT\xa0\xa09.02\xa0(60663832)\xa0\xa0\xa0AAT\xa018.26\xa0(122860591)\xa0\xa0AGT\xa0\xa09.08\xa0(61127436)\xa0\n',
                              'ATC\xa024.60\xa0(165520671)\xa0\xa0ACC\xa022.88\xa0(153963097)\xa0\xa0AAC\xa021.47\xa0(144468591)\xa0\xa0AGC\xa015.89\xa0(106898991)\n',
                              'ATA\xa0\xa04.88\xa0(32821498)\xa0\xa0\xa0ACA\xa0\xa07.63\xa0(51354438)\xa0\xa0\xa0AAA\xa033.94\xa0(228397925)\xa0\xa0AGA\xa0\xa02.43\xa0(16341387)\xa0\n',
                              'ATG\xa027.59\xa0(185678714)\xa0\xa0ACG\xa014.47\xa0(97383052)\xa0\xa0\xa0AAG\xa010.70\xa0(71980179)\xa0\xa0\xa0AGG\xa0\xa01.48\xa0(9931996)\xa0\xa0\n', 
                              '\n', 
                              'GTT\xa018.39\xa0(123762454)\xa0\xa0GCT\xa015.54\xa0(104562996)\xa0\xa0GAT\xa032.43\xa0(218213146)\xa0\xa0GGT\xa024.45\xa0(164526277)\n', 
                              'GTC\xa015.07\xa0(101415663)\xa0\xa0GCC\xa025.45\xa0(171276274)\xa0\xa0GAC\xa019.14\xa0(128814139)\xa0\xa0GGC\xa028.65\xa0(192822620)\n',
                              'GTA\xa010.97\xa0(73847538)\xa0\xa0\xa0GCA\xa020.61\xa0(138682074)\xa0\xa0GAA\xa039.55\xa0(266144484)\xa0\xa0GGA\xa0\xa08.44\xa0(56771393)\xa0\n',
                              'GTG\xa025.90\xa0(174290187)\xa0\xa0GCG\xa032.79\xa0(220665389)\xa0\xa0GAG\xa018.24\xa0(122743123)\xa0\xa0GGG\xa011.29\xa0(75943843)\xa0'],


        'Homo_sapien' : ['Table contains 159069 CDSs (99738040 codons)\n',
                         '\n',
                         'TTT\xa017.06\xa0(1701077)\xa0\xa0TCT\xa016.58\xa0(1653836)\xa0\xa0TAT\xa012.04\xa0(1200411)\xa0\xa0TGT\xa010.54\xa0(1050785)\n',
                         'TTC\xa017.87\xa0(1782473)\xa0\xa0TCC\xa017.44\xa0(1739446)\xa0\xa0TAC\xa013.70\xa0(1366156)\xa0\xa0TGC\xa011.15\xa0(1112221)\n',
                         'TTA\xa0\xa08.55\xa0(853143)\xa0\xa0\xa0TCA\xa013.89\xa0(1385383)\xa0\xa0TAA\xa0\xa00.46\xa0(45381)\xa0\xa0\xa0\xa0TGA\xa0\xa00.83\xa0(83155)\xa0\xa0\n',
                         'TTG\xa013.30\xa0(1326294)\xa0\xa0TCG\xa0\xa04.18\xa0(417288)\xa0\xa0\xa0TAG\xa0\xa00.36\xa0(36112)\xa0\xa0\xa0\xa0TGG\xa011.77\xa0(1173774)\n', 
                         '\n', 
                         'CTT\xa013.95\xa0(1390970)\xa0\xa0CCT\xa018.88\xa0(1882918)\xa0\xa0CAT\xa011.74\xa0(1170833)\xa0\xa0CGT\xa0\xa04.54\xa0(452912)\xa0\n',
                         'CTC\xa018.06\xa0(1801531)\xa0\xa0CCC\xa019.19\xa0(1913483)\xa0\xa0CAC\xa014.76\xa0(1472070)\xa0\xa0CGC\xa0\xa09.06\xa0(904107)\xa0\n',
                         'CTA\xa0\xa07.39\xa0(736716)\xa0\xa0\xa0CCA\xa018.45\xa0(1839881)\xa0\xa0CAA\xa013.83\xa0(1379401)\xa0\xa0CGA\xa0\xa06.36\xa0(634193)\xa0\n',
                         'CTG\xa036.75\xa0(3665034)\xa0\xa0CCG\xa0\xa06.36\xa0(633907)\xa0\xa0\xa0CAG\xa035.30\xa0(3520920)\xa0\xa0CGG\xa010.88\xa0(1085602)\n',
                         '\n', 
                         'ATT\xa016.36\xa0(1631789)\xa0\xa0ACT\xa014.12\xa0(1408630)\xa0\xa0AAT\xa018.16\xa0(1810933)\xa0\xa0AGT\xa013.72\xa0(1368732)\n',
                         'ATC\xa018.97\xa0(1891991)\xa0\xa0ACC\xa017.95\xa0(1790317)\xa0\xa0AAC\xa018.36\xa0(1831427)\xa0\xa0AGC\xa019.74\xa0(1969101)\n', 
                         'ATA\xa0\xa07.98\xa0(796217)\xa0\xa0\xa0ACA\xa016.33\xa0(1629094)\xa0\xa0AAA\xa027.15\xa0(2708089)\xa0\xa0AGA\xa013.09\xa0(1305261)\n', 
                         'ATG\xa021.40\xa0(2134650)\xa0\xa0ACG\xa0\xa05.72\xa0(570644)\xa0\xa0\xa0AAG\xa031.89\xa0(3180910)\xa0\xa0AGG\xa012.15\xa0(1212181)\n',
                         '\n', 
                         'GTT\xa011.59\xa0(1155684)\xa0\xa0GCT\xa018.77\xa0(1872140)\xa0\xa0GAT\xa023.68\xa0(2361597)\xa0\xa0GGT\xa010.75\xa0(1072606)\n',
                         'GTC\xa013.58\xa0(1354925)\xa0\xa0GCC\xa026.18\xa0(2611615)\xa0\xa0GAC\xa024.49\xa0(2442300)\xa0\xa0GGC\xa020.23\xa0(2017514)\n',
                         'GTA\xa0\xa07.56\xa0(753779)\xa0\xa0\xa0GCA\xa016.89\xa0(1684588)\xa0\xa0GAA\xa033.04\xa0(3294994)\xa0\xa0GGA\xa017.02\xa0(1697657)\n',
                         'GTG\xa026.24\xa0(2616674)\xa0\xa0GCG\xa0\xa06.26\xa0(624561)\xa0\xa0\xa0GAG\xa039.88\xa0(3977521)\xa0\xa0GGG\xa015.53\xa0(1548506)'],


        'Mus_musculus' : ['Table contains 78862 CDSs (53348424 codons)\n',
              '\n',
              'TTT\xa015.94\xa0(850518)\xa0\xa0\xa0TCT\xa017.39\xa0(927915)\xa0\xa0\xa0TAT\xa011.15\xa0(594603)\xa0\xa0\xa0TGT\xa010.68\xa0(569800)\xa0\n',
              'TTC\xa018.81\xa0(1003402)\xa0\xa0TCC\xa018.32\xa0(977319)\xa0\xa0\xa0TAC\xa014.42\xa0(769179)\xa0\xa0\xa0TGC\xa010.95\xa0(584071)\xa0\n',
              'TTA\xa0\xa07.29\xa0(388732)\xa0\xa0\xa0TCA\xa013.31\xa0(710152)\xa0\xa0\xa0TAA\xa0\xa00.39\xa0(20957)\xa0\xa0\xa0\xa0TGA\xa0\xa00.76\xa0(40417)\xa0\xa0\n', 
              'TTG\xa013.27\xa0(708046)\xa0\xa0\xa0TCG\xa0\xa04.29\xa0(228719)\xa0\xa0\xa0TAG\xa0\xa00.35\xa0(18804)\xa0\xa0\xa0\xa0TGG\xa011.44\xa0(610420)\xa0\n',
              '\n',
              'CTT\xa013.45\xa0(717284)\xa0\xa0\xa0CCT\xa020.06\xa0(1070067)\xa0\xa0CAT\xa011.23\xa0(599053)\xa0\xa0\xa0CGT\xa0\xa04.64\xa0(247431)\xa0\n',
              'CTC\xa018.83\xa0(1004744)\xa0\xa0CCC\xa018.34\xa0(978608)\xa0\xa0\xa0CAC\xa015.23\xa0(812332)\xa0\xa0\xa0CGC\xa0\xa08.60\xa0(459061)\xa0\n',
              'CTA\xa0\xa08.04\xa0(428925)\xa0\xa0\xa0CCA\xa019.05\xa0(1016527)\xa0\xa0CAA\xa013.11\xa0(699160)\xa0\xa0\xa0CGA\xa0\xa06.90\xa0(368006)\xa0\n',
              'CTG\xa037.31\xa0(1990552)\xa0\xa0CCG\xa0\xa06.11\xa0(326210)\xa0\xa0\xa0CAG\xa036.71\xa0(1958336)\xa0\xa0CGG\xa010.46\xa0(558268)\xa0\n',
              '\n',
              'ATT\xa014.66\xa0(782173)\xa0\xa0\xa0ACT\xa013.92\xa0(742744)\xa0\xa0\xa0AAT\xa015.90\xa0(848342)\xa0\xa0\xa0AGT\xa014.11\xa0(752647)\xa0\n',
              'ATC\xa020.33\xa0(1084325)\xa0\xa0ACC\xa018.16\xa0(968779)\xa0\xa0\xa0AAC\xa019.75\xa0(1053455)\xa0\xa0AGC\xa020.77\xa0(1107848)\n', 
              'ATA\xa0\xa07.28\xa0(388619)\xa0\xa0\xa0ACA\xa016.67\xa0(889155)\xa0\xa0\xa0AAA\xa023.84\xa0(1271932)\xa0\xa0AGA\xa013.03\xa0(695038)\xa0\n',
              'ATG\xa021.70\xa0(1157638)\xa0\xa0ACG\xa0\xa05.73\xa0(305495)\xa0\xa0\xa0AAG\xa033.95\xa0(1811037)\xa0\xa0AGG\xa012.94\xa0(690288)\xa0\n',
              '\n',
              'GTT\xa010.81\xa0(576778)\xa0\xa0\xa0GCT\xa020.19\xa0(1077009)\xa0\xa0GAT\xa022.33\xa0(1191101)\xa0\xa0GGT\xa011.09\xa0(591536)\xa0\n',
              'GTC\xa014.52\xa0(774776)\xa0\xa0\xa0GCC\xa025.16\xa0(1342410)\xa0\xa0GAC\xa026.30\xa0(1403166)\xa0\xa0GGC\xa019.81\xa0(1056977)\n',
              'GTA\xa0\xa07.48\xa0(399065)\xa0\xa0\xa0GCA\xa016.80\xa0(896415)\xa0\xa0\xa0GAA\xa030.33\xa0(1618320)\xa0\xa0GGA\xa016.77\xa0(894696)\xa0\n',
              'GTG\xa026.58\xa0(1417892)\xa0\xa0GCG\xa0\xa05.86\xa0(312882)\xa0\xa0\xa0GAG\xa041.48\xa0(2213104)\xa0\xa0GGG\xa014.91\xa0(795164)\xa0'],


        'Caenorhabditis_elegans' : ['Table contains 28093 CDSs (12891394 codons)\n',
                                    '\n',
                                    'TTT\xa021.42\xa0(276101)\xa0\xa0TCT\xa016.96\xa0(218675)\xa0\xa0TAT\xa017.11\xa0(220515)\xa0\xa0TGT\xa010.99\xa0(141670)\n',
                                    'TTC\xa023.05\xa0(297187)\xa0\xa0TCC\xa010.61\xa0(136826)\xa0\xa0TAC\xa013.48\xa0(173738)\xa0\xa0TGC\xa0\xa08.76\xa0(112920)\n',
                                    'TTA\xa0\xa09.31\xa0(120042)\xa0\xa0TCA\xa021.14\xa0(272536)\xa0\xa0TAA\xa0\xa01.07\xa0(13817)\xa0\xa0\xa0TGA\xa0\xa00.72\xa0(9285)\xa0\xa0\n',
                                    'TTG\xa019.58\xa0(252430)\xa0\xa0TCG\xa012.81\xa0(165145)\xa0\xa0TAG\xa0\xa00.39\xa0(4993)\xa0\xa0\xa0\xa0TGG\xa010.68\xa0(137704)\n',
                                    '\n',
                                    'CTT\xa020.92\xa0(269646)\xa0\xa0CCT\xa0\xa09.16\xa0(118092)\xa0\xa0CAT\xa014.30\xa0(184296)\xa0\xa0CGT\xa011.38\xa0(146661)\n', 
                                    'CTC\xa014.54\xa0(187429)\xa0\xa0CCC\xa0\xa04.46\xa0(57433)\xa0\xa0\xa0CAC\xa0\xa09.07\xa0(116940)\xa0\xa0CGC\xa0\xa05.03\xa0(64816)\xa0\n',
                                    'CTA\xa0\xa07.69\xa0(99108)\xa0\xa0\xa0CCA\xa027.25\xa0(351282)\xa0\xa0CAA\xa027.88\xa0(359465)\xa0\xa0CGA\xa012.50\xa0(161206)\n',
                                    'CTG\xa012.05\xa0(155316)\xa0\xa0CCG\xa010.28\xa0(132512)\xa0\xa0CAG\xa014.88\xa0(191850)\xa0\xa0CGG\xa0\xa04.84\xa0(62382)\xa0\n',
                                    '\n',
                                    'ATT\xa031.51\xa0(406147)\xa0\xa0ACT\xa019.30\xa0(248798)\xa0\xa0AAT\xa030.06\xa0(387496)\xa0\xa0AGT\xa012.28\xa0(158278)\n', 
                                    'ATC\xa018.50\xa0(238544)\xa0\xa0ACC\xa010.32\xa0(132986)\xa0\xa0AAC\xa017.93\xa0(231169)\xa0\xa0AGC\xa0\xa08.32\xa0(107288)\n',
                                    'ATA\xa0\xa09.03\xa0(116402)\xa0\xa0ACA\xa020.55\xa0(264946)\xa0\xa0AAA\xa036.40\xa0(469296)\xa0\xa0AGA\xa015.25\xa0(196646)\n', 
                                    'ATG\xa026.09\xa0(336309)\xa0\xa0ACG\xa0\xa09.22\xa0(118904)\xa0\xa0AAG\xa025.58\xa0(329818)\xa0\xa0AGG\xa0\xa03.80\xa0(48928)\xa0\n',
                                    '\n',
                                    'GTT\xa024.14\xa0(311209)\xa0\xa0GCT\xa022.94\xa0(295696)\xa0\xa0GAT\xa036.73\xa0(473446)\xa0\xa0GGT\xa011.12\xa0(143295)\n',
                                    'GTC\xa013.58\xa0(175030)\xa0\xa0GCC\xa012.74\xa0(164252)\xa0\xa0GAC\xa017.27\xa0(222661)\xa0\xa0GGC\xa0\xa06.78\xa0(87341)\xa0\n', 
                                    'GTA\xa0\xa09.78\xa0(126083)\xa0\xa0GCA\xa020.35\xa0(262322)\xa0\xa0GAA\xa041.61\xa0(536431)\xa0\xa0GGA\xa032.03\xa0(412865)\n',
                                    'GTG\xa014.53\xa0(187301)\xa0\xa0GCG\xa0\xa08.56\xa0(110317)\xa0\xa0GAG\xa025.10\xa0(323529)\xa0\xa0GGG\xa0\xa04.32\xa0(55643)\xa0'],

    
    
    'Saccharomyces_cerevisiae' : ['Table contains 5888 CDSs (2917768 codons)\n',
                                      '\n',
                                      'TTT\xa026.18\xa0(76390)\xa0\xa0TCT\xa023.35\xa0(68138)\xa0\xa0TAT\xa019.05\xa0(55593)\xa0\xa0\xa0TGT\xa0\xa07.82\xa0(22826)\n',
                                      'TTC\xa017.88\xa0(52158)\xa0\xa0TCC\xa014.07\xa0(41039)\xa0\xa0TAC\xa014.60\xa0(42591)\xa0\xa0\xa0TGC\xa0\xa04.75\xa0(13847)\n', 
                                      'TTA\xa026.33\xa0(76815)\xa0\xa0TCA\xa019.05\xa0(55589)\xa0\xa0TAA\xa0\xa00.95\xa0(2786)\xa0\xa0\xa0\xa0TGA\xa0\xa00.60\xa0(1750)\xa0\n',
                                      'TTG\xa026.50\xa0(77319)\xa0\xa0TCG\xa0\xa08.71\xa0(25412)\xa0\xa0TAG\xa0\xa00.46\xa0(1352)\xa0\xa0\xa0\xa0TGG\xa010.35\xa0(30212)\n', 
                                      '\n',
                                      'CTT\xa012.27\xa0(35793)\xa0\xa0CCT\xa013.57\xa0(39595)\xa0\xa0CAT\xa013.89\xa0(40523)\xa0\xa0\xa0CGT\xa0\xa06.26\xa0(18269)\n',
                                      'CTC\xa0\xa05.52\xa0(16115)\xa0\xa0CCC\xa0\xa06.91\xa0(20160)\xa0\xa0CAC\xa0\xa07.74\xa0(22597)\xa0\xa0\xa0CGC\xa0\xa02.63\xa0(7661)\xa0\n',
                                      'CTA\xa013.52\xa0(39441)\xa0\xa0CCA\xa017.81\xa0(51970)\xa0\xa0CAA\xa027.10\xa0(79082)\xa0\xa0\xa0CGA\xa0\xa03.10\xa0(9053)\xa0\n',
                                      'CTG\xa010.65\xa0(31077)\xa0\xa0CCG\xa0\xa05.42\xa0(15808)\xa0\xa0CAG\xa012.42\xa0(36228)\xa0\xa0\xa0CGG\xa0\xa01.82\xa0(5310)\xa0\n',
                                      '\n',
                                      'ATT\xa030.10\xa0(87820)\xa0\xa0ACT\xa020.24\xa0(59052)\xa0\xa0AAT\xa036.61\xa0(106809)\xa0\xa0AGT\xa014.60\xa0(42600)\n',
                                      'ATC\xa016.99\xa0(49563)\xa0\xa0ACC\xa012.48\xa0(36428)\xa0\xa0AAC\xa024.80\xa0(72354)\xa0\xa0\xa0AGC\xa0\xa09.96\xa0(29053)\n',
                                      'ATA\xa018.29\xa0(53372)\xa0\xa0ACA\xa018.18\xa0(53050)\xa0\xa0AAA\xa042.83\xa0(124954)\xa0\xa0AGA\xa021.05\xa0(61405)\n', 
                                      'ATG\xa020.68\xa0(60332)\xa0\xa0ACG\xa0\xa08.15\xa0(23785)\xa0\xa0AAG\xa030.52\xa0(89048)\xa0\xa0\xa0AGG\xa0\xa09.45\xa0(27577)\n', 
                                      '\n',
                                      'GTT\xa021.47\xa0(62642)\xa0\xa0GCT\xa020.28\xa0(59162)\xa0\xa0GAT\xa038.09\xa0(111152)\xa0\xa0GGT\xa022.59\xa0(65916)\n', 
                                      'GTC\xa011.23\xa0(32773)\xa0\xa0GCC\xa012.14\xa0(35414)\xa0\xa0GAC\xa020.39\xa0(59500)\xa0\xa0\xa0GGC\xa0\xa09.78\xa0(28550)\n', 
                                      'GTA\xa012.07\xa0(35204)\xa0\xa0GCA\xa016.26\xa0(47457)\xa0\xa0GAA\xa045.81\xa0(133662)\xa0\xa0GGA\xa011.19\xa0(32639)\n',
                                      'GTG\xa010.72\xa0(31271)\xa0\xa0GCG\xa0\xa06.17\xa0(17996)\xa0\xa0GAG\xa019.55\xa0(57048)\xa0\xa0\xa0GGG\xa0\xa06.06\xa0(17681)']
}


speciesList = ['Escherichia coli', 'Caenorhabditis elegans', 'Mus musculus', 'Homo sapien', 'Saccharomyces cerevisiae']
speciesList2 = ['escherichia coli', 'caenorhabditis elegans', 'mus musculus', 'homo sapien', 'saccharomyces cerevisiae', '1', '2', '3', '4','5']
#Allows numbers or names to be entered for ease of use
inputDict = {'escherichia coli':'Escherichia_coli', '1': 'Escherichia_coli', 'caenorhabditis elegans': 'Caenorhabditis_elegans',
             '2':'Caenorhabditis_elegans', 'mus musculus':'Mus_musculus', '3': 'Mus_musculus','homo sapien':'Homo_sapien',
             '4':'Homo_sapien', 'saccharomyces cerevisiae':'Saccharomyces_cerevisiae', '5':'Saccharomyces_cerevisiae'}
#does MinMax calculation
def calculateMinMax(sequence, aaFreqDict, freqDict, mapDict, windowSize):
    freqDict = freqDict
    aaFreqDict = aaFreqDict
    windowSize = windowSize
    mapDict = mapDict
    codonSeq = sequence
    minMaxValues = []
    
    
    for i in range(int(windowSize/2)):
        minMaxValues.append(0)
    
    #Using the specified sliding window size (windowSize/2 - 1 on either side of the central codon), min/max is calculated
    for i in range(len(codonSeq)-windowSize+1):
        window = codonSeq[i:i+windowSize] #list of the codons in the current window

        Actual = 0.0     #average of the actual codon frequencies
        Max = 0.0        #average of the min codon frequencies
        Min = 0.0        #average of the max codon frequencies
        Avg = 0.0        #average of the averages of all the frequencies associated with each amino acid

        #Sum the frequencies
        for codon in window:
            frequencies = aaFreqDict[mapDict[codon]] #list of all frequencies associated with the amino acid this codon encodes

            Actual += freqDict[codon]
            Max += max(frequencies)
            Min += min(frequencies)
            Avg += sum(frequencies)/len(frequencies)

        #Divide by the window size to get the averages
        Actual = Actual/windowSize
        Max = Max/windowSize
        Min = Min/windowSize
        Avg = Avg/windowSize

        percentMax = ((Actual-Avg)/(Max-Avg))*100
        percentMin = ((Avg-Actual)/(Avg-Min))*100

        if(percentMax >= 0):
            minMaxValues.append(round(percentMax,2))
        else:
            minMaxValues.append(round(-percentMin,2))

    #fills in values for codons where window size makes min/max unable to be calculated
    for i in range(int(windowSize/2)):
        minMaxValues.append(0)

    return minMaxValues

def sign(a):
	return (1+a)/(1+abs(a))

def fitness(minMaxValues,aaFreqDict, freqDict, mapDict, windowSize):
	def g(solution):
		nowvals = calculateMinMax(solution.dna(),aaFreqDict, freqDict, mapDict, windowSize)
		sum_error = 0
		maxval = 0
		for k in range(0,len(nowvals)):
			val = abs(nowvals[k] - minMaxValues[k])
			if k > 1 and not sign(minMaxValues[k] - minMaxValues[k-1]) == sign(nowvals[k]-nowvals[k-1]):
				sum_error += 100*val
			else:
				sum_error += val
			if val > maxval:
				maxval = val
		return sum_error + maxval
	return g
def mmin(a):
	p = (2**32,0)
	for y in a:
		if y[0] < p[0]:
			p = y
	return p
def worst_genes(minMaxValues, aaFreqDict, freqDict, MapDict, WindowSize):
	def k(solution):
		nowvals = calculateMinMax(solution.dna(),aaFreqDict, freqDict, mapDict,windowSize)
		worst = [(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0)]
		for k in range(0,len(nowvals)):
			sm = mmin(worst)
			val = abs(nowvals[k] - minMaxValues[k]) 
			if val > sm[0]:
				worst.remove(sm)
				worst.append((val,k))
		mworst = set()
		for k in worst:
			mval = k[1]
			for y in range(int(mval-(windowSize/2)),int(mval+(windowSize/2)+1)):
				mworst.add(y)
		return list(mworst)
	return k




#print("You may either use one of the following preloaded codon frequency tables")
#for i in range(int(len(speciesList))):
#    print(i+1, speciesList[i])
usageFile = "1" #input("Input species for usage file (name or number above): ")
usageFile = usageFile.lower()
while usageFile not in speciesList2:
    print("")
    print("Unrecognized input, please enter either species name or index number: ")
    for i in range(int(len(speciesList2)/2)):
        print(i+1, speciesList[i])
    usageFile = input("Input species for usage file (name or corresponding number): ")
    usageFile = usageFile.lower()
frequenciesFile = speciesDict[inputDict[usageFile]]
    
#Sequence Input Handling    
sequence = "ATGCTTAAACAGGTAGAAATTTTCACCGATGGTTCGTGTCTGGGCAATCCAGGACCTGGGGGTTACGGCGCTATTTTACGCTATCGCGGACGCGAGAAAACCTTTAGCGCTGGCTACACCCGCACCACCAACAACCGTATGGAGTTGATGGCCGCTATTGTCGCGCTGGAGGCGTTAAAAGAACATTGCGAAGTCATTTTGAGTACCGACAGCCAGTATGTCCGCCAGGGTATCACCCAGTGGATCCATAACTGGAAAAAACGTGGCTGGAAAACCGCAGACAAAAAACCAGTAAAAAATGTCGATCTCTGGCAACGTCTTGATGCTGCATTGGGGCAGCATCAAATCAAATGGGAATGGGTTAAAGGCCATGCCGGACACCCGGAAAACGAACGCTGTGATGAACTGGCTCGTGCCGCGGCGATGAATCCCACACTGGAAGATACAGGCTACCAAGTTGAAGTT"#input('Nucleotide Sequence (start to stop codon): ')
sequence = sequence.upper()
if sequence[0:3] != "ATG":
    response = input("Sequence does not begin with ATG, would you re-enter the sequence (yes or no): ")
    response = response.lower()
    while response == "yes":
        sequence = input("Please reenter the codon sequence: ")
        if sequence[0:3] != "ATG":
            response = input("Sequence does not begin with ATG, would you like to re-enter the sequence (yes or no): ")
        else:
            response = "no"
while len(sequence)%3 != 0:
    sequence = input("The entered nucleotide sequence is not divisible by three, please reenter the sequence:")
    sequence = sequence.upper()

print("Calculating %MinMax...")

#Data Cleaning
for line in frequenciesFile:
    line = line.split()
    i=0
    if len(line)>11:  
        while i < len(line):
            freqDict[line[i]] = float(line[i+1])
            aaMapDict[mapDict[line[i]]] = []
            aaFreqDict[mapDict[line[i]]] = []
            i+=3

for line in frequenciesFile:
    line = line.split()
    i = 0
    if len(line)>11:
        while i<len(line):
            aaFreqDict[mapDict[line[i]]].append(float(line[i+1]))
            aaMapDict[mapDict[line[i]]].append(str(line[i] + " " + line[i+1]))
            i+=3
    
#For a given input fasta sequence, break into codons and corresponding amino acids
codonSeq = []

extras = ""
for line in sequence:
    line = line.rstrip()
    string = str(extras) + str(line)
    i=0
    j=3
    while j<=len(string):
        codonSeq.append(string[i:j])
        i+=3
        j+=3
    extras = str(string[i:])

aaSeq = []
for codon in codonSeq:
    aaSeq.append(mapDict[codon])
     
minMaxValues = calculateMinMax(codonSeq, aaFreqDict, freqDict, mapDict, windowSize)




################################################################ MPI #########################################################################

if rank > 0:
    while True:
        data = comm.recv(source=0,tag=1)
        #print(rank," Got data")
        child = mpiAlgo.multiproc_vertex(data)
        #print(child)
        comm.send(child,dest=0,tag=1)
        #req.wait()

################################################################ MPI ##########################################################################



#########################################     5    ###################################
usageFile2 = "5" #input("Input species for harmonization: ")
usageFile2 = usageFile2.lower()
while usageFile2 not in speciesList2:
    print("")
    print("Unrecognized input, please enter either species name or index number: ")
    for i in range(int(len(speciesList2)/2)):
        print(i+1, speciesList[i])
    usageFile2 = input("Input species for usage file (name or corresponding number): ")
    usageFile2 = usageFile2.lower()
frequenciesFile2 = speciesDict[inputDict[usageFile2]]

freqDict2 = dict()
aaMapDict2 = dict()#dictionary from amino acid to list of codons with frequencies for it (for RRTs)
aaFreqDict2 = dict()
for line in frequenciesFile2:
    line = line.split()
    i=0
    if len(line)>11:
        while i < len(line):
            freqDict2[line[i]] = float(line[i+1])
            aaMapDict2[mapDict[line[i]]] = []
            aaFreqDict2[mapDict[line[i]]] = []
            i+=3

for line in frequenciesFile2:
    line = line.split()
    i = 0
    if len(line)>11:
        while i<len(line):
            aaFreqDict2[mapDict[line[i]]].append(float(line[i+1]))
            aaMapDict2[mapDict[line[i]]].append(str(line[i] + " " + line[i+1]))
            i+=3


mfitness = fitness(minMaxValues, aaFreqDict2, freqDict2, mapDict, windowSize)
mworst = worst_genes(minMaxValues, aaFreqDict2, freqDict2, mapDict, windowSize)
start = [random.choice(mpiAlgo.rev[mpiAlgo.mapDict[x]]) for x in codonSeq]

a = []
b = []
for k in range(10):
    print("Num 5, loop %i"%k)
    best = mpiAlgo.graph_run(mpiAlgo.solution(start),mfitness,mworst,10,100,50)
    a.append(best[2])
    b.append(best[1])
    #print(best)
print("Num 5 Time Mean: %f stdev %f"%(mean(a),stdev(a)))
print("Num 5 Score Mean: %f stdev %f"%(mean(b),stdev(b)))

##print(best)

orig = open('orig.dat','w+')
for k in range(0,len(minMaxValues)):
	orig.write("\n"+str(k)+" "+str(minMaxValues[k]))
orig.close()
new = open('new-5.dat','w+')
newScores = calculateMinMax(best[0],aaFreqDict2, freqDict2, mapDict, windowSize)
for k in range(0,len(newScores)):
	new.write("\n"+str(k)+" "+str(newScores[k]))
new.close()


#########################################     4    ###################################
usageFile2 = "4" #input("Input species for harmonization: ")
usageFile2 = usageFile2.lower()
while usageFile2 not in speciesList2:
    print("")
    print("Unrecognized input, please enter either species name or index number: ")
    for i in range(int(len(speciesList2)/2)):
        print(i+1, speciesList[i])
    usageFile2 = input("Input species for usage file (name or corresponding number): ")
    usageFile2 = usageFile2.lower()
frequenciesFile2 = speciesDict[inputDict[usageFile2]]

freqDict2 = dict()
aaMapDict2 = dict()#dictionary from amino acid to list of codons with frequencies for it (for RRTs)
aaFreqDict2 = dict()
for line in frequenciesFile2:
    line = line.split()
    i=0
    if len(line)>11:
        while i < len(line):
            freqDict2[line[i]] = float(line[i+1])
            aaMapDict2[mapDict[line[i]]] = []
            aaFreqDict2[mapDict[line[i]]] = []
            i+=3

for line in frequenciesFile2:
    line = line.split()
    i = 0
    if len(line)>11:
        while i<len(line):
            aaFreqDict2[mapDict[line[i]]].append(float(line[i+1]))
            aaMapDict2[mapDict[line[i]]].append(str(line[i] + " " + line[i+1]))
            i+=3

mfitness = fitness(minMaxValues, aaFreqDict2, freqDict2, mapDict, windowSize)
mworst = worst_genes(minMaxValues, aaFreqDict2, freqDict2, mapDict, windowSize)
start = [random.choice(mpiAlgo.rev[mpiAlgo.mapDict[x]]) for x in codonSeq]

a = []
b = []
for k in range(10):
    print("Num 4, loop %i"%k)
    best = mpiAlgo.graph_run(mpiAlgo.solution(start),mfitness,mworst,10,100,50)
    a.append(best[2])
    b.append(best[1])
    #print(best)
print("Num 4 Time Mean: %f stdev %f"%(mean(a),stdev(a)))
print("Num 4 Score Mean: %f stdev %f"%(mean(b),stdev(b)))




orig = open('orig.dat','w+')
for k in range(0,len(minMaxValues)):
	orig.write("\n"+str(k)+" "+str(minMaxValues[k]))
orig.close()
new = open('new-4.dat','w+')
newScores = calculateMinMax(best[0],aaFreqDict2, freqDict2, mapDict, windowSize)
for k in range(0,len(newScores)):
	new.write("\n"+str(k)+" "+str(newScores[k]))
new.close()

#########################################     3    ###################################
usageFile2 = "3" #input("Input species for harmonization: ")
usageFile2 = usageFile2.lower()
while usageFile2 not in speciesList2:
    print("")
    print("Unrecognized input, please enter either species name or index number: ")
    for i in range(int(len(speciesList2)/2)):
        print(i+1, speciesList[i])
    usageFile2 = input("Input species for usage file (name or corresponding number): ")
    usageFile2 = usageFile2.lower()
frequenciesFile2 = speciesDict[inputDict[usageFile2]]

freqDict2 = dict()
aaMapDict2 = dict()#dictionary from amino acid to list of codons with frequencies for it (for RRTs)
aaFreqDict2 = dict()
for line in frequenciesFile2:
    line = line.split()
    i=0
    if len(line)>11:
        while i < len(line):
            freqDict2[line[i]] = float(line[i+1])
            aaMapDict2[mapDict[line[i]]] = []
            aaFreqDict2[mapDict[line[i]]] = []
            i+=3

for line in frequenciesFile2:
    line = line.split()
    i = 0
    if len(line)>11:
        while i<len(line):
            aaFreqDict2[mapDict[line[i]]].append(float(line[i+1]))
            aaMapDict2[mapDict[line[i]]].append(str(line[i] + " " + line[i+1]))
            i+=3

mfitness = fitness(minMaxValues, aaFreqDict2, freqDict2, mapDict, windowSize)
mworst = worst_genes(minMaxValues, aaFreqDict2, freqDict2, mapDict, windowSize)
start = [random.choice(mpiAlgo.rev[mpiAlgo.mapDict[x]]) for x in codonSeq]

a = []
b = []
for k in range(10):
    print("Num 3, loop %i"%k)
    best = mpiAlgo.graph_run(mpiAlgo.solution(start),mfitness,mworst,10,100,50)
    a.append(best[2])
    b.append(best[1])
    #print(best)
print("Num 3 Time Mean: %f stdev %f"%(mean(a),stdev(a)))
print("Num 3 Score Mean: %f stdev %f"%(mean(b),stdev(b)))


orig = open('orig.dat','w+')
for k in range(0,len(minMaxValues)):
	orig.write("\n"+str(k)+" "+str(minMaxValues[k]))
orig.close()
new = open('new-3.dat','w+')
newScores = calculateMinMax(best[0],aaFreqDict2, freqDict2, mapDict, windowSize)
for k in range(0,len(newScores)):
	new.write("\n"+str(k)+" "+str(newScores[k]))
new.close()

#########################################     2    ###################################
usageFile2 = "2" #input("Input species for harmonization: ")
usageFile2 = usageFile2.lower()
while usageFile2 not in speciesList2:
    print("")
    print("Unrecognized input, please enter either species name or index number: ")
    for i in range(int(len(speciesList2)/2)):
        print(i+1, speciesList[i])
    usageFile2 = input("Input species for usage file (name or corresponding number): ")
    usageFile2 = usageFile2.lower()
frequenciesFile2 = speciesDict[inputDict[usageFile2]]

freqDict2 = dict()
aaMapDict2 = dict()#dictionary from amino acid to list of codons with frequencies for it (for RRTs)
aaFreqDict2 = dict()
for line in frequenciesFile2:
    line = line.split()
    i=0
    if len(line)>11:
        while i < len(line):
            freqDict2[line[i]] = float(line[i+1])
            aaMapDict2[mapDict[line[i]]] = []
            aaFreqDict2[mapDict[line[i]]] = []
            i+=3

for line in frequenciesFile2:
    line = line.split()
    i = 0
    if len(line)>11:
        while i<len(line):
            aaFreqDict2[mapDict[line[i]]].append(float(line[i+1]))
            aaMapDict2[mapDict[line[i]]].append(str(line[i] + " " + line[i+1]))
            i+=3

mfitness = fitness(minMaxValues, aaFreqDict2, freqDict2, mapDict, windowSize)
mworst = worst_genes(minMaxValues, aaFreqDict2, freqDict2, mapDict, windowSize)
start = [random.choice(mpiAlgo.rev[mpiAlgo.mapDict[x]]) for x in codonSeq]


a = []
b = []
for k in range(10):
    print("Num 2, loop %i"%k)
    best = mpiAlgo.graph_run(mpiAlgo.solution(start),mfitness,mworst,10,100,50)
    a.append(best[2])
    b.append(best[1])
    #print(best)
print("Num 2 Time Mean: %f stdev %f"%(mean(a),stdev(a)))
print("Num 2 Score Mean: %f stdev %f"%(mean(b),stdev(b)))

orig = open('orig.dat','w+')
for k in range(0,len(minMaxValues)):
	orig.write("\n"+str(k)+" "+str(minMaxValues[k]))
orig.close()
new = open('new-2.dat','w+')
newScores = calculateMinMax(best[0],aaFreqDict2, freqDict2, mapDict, windowSize)
for k in range(0,len(newScores)):
	new.write("\n"+str(k)+" "+str(newScores[k]))
new.close()

#########################################     1    ###################################
usageFile2 = "1" #input("Input species for harmonization: ")
usageFile2 = usageFile2.lower()
while usageFile2 not in speciesList2:
    print("")
    print("Unrecognized input, please enter either species name or index number: ")
    for i in range(int(len(speciesList2)/2)):
        print(i+1, speciesList[i])
    usageFile2 = input("Input species for usage file (name or corresponding number): ")
    usageFile2 = usageFile2.lower()
frequenciesFile2 = speciesDict[inputDict[usageFile2]]

freqDict2 = dict()
aaMapDict2 = dict()#dictionary from amino acid to list of codons with frequencies for it (for RRTs)
aaFreqDict2 = dict()
for line in frequenciesFile2:
    line = line.split()
    i=0
    if len(line)>11:
        while i < len(line):
            freqDict2[line[i]] = float(line[i+1])
            aaMapDict2[mapDict[line[i]]] = []
            aaFreqDict2[mapDict[line[i]]] = []
            i+=3

for line in frequenciesFile2:
    line = line.split()
    i = 0
    if len(line)>11:
        while i<len(line):
            aaFreqDict2[mapDict[line[i]]].append(float(line[i+1]))
            aaMapDict2[mapDict[line[i]]].append(str(line[i] + " " + line[i+1]))
            i+=3

mfitness = fitness(minMaxValues, aaFreqDict2, freqDict2, mapDict, windowSize)
mworst = worst_genes(minMaxValues, aaFreqDict2, freqDict2, mapDict, windowSize)
start = [random.choice(mpiAlgo.rev[mpiAlgo.mapDict[x]]) for x in codonSeq]


a = []
b = []
for k in range(10):
    print("Num 1, loop %i"%k)
    best = mpiAlgo.graph_run(mpiAlgo.solution(start),mfitness,mworst,10,100,50)
    a.append(best[2])
    b.append(best[1])
    #print(best)
print("Num 1 Time Mean: %f stdev %f"%(mean(a),stdev(a)))
print("Num 1 Score Mean: %f stdev %f"%(mean(b),stdev(b)))


orig = open('orig.dat','w+')
for k in range(0,len(minMaxValues)):
	orig.write("\n"+str(k)+" "+str(minMaxValues[k]))
orig.close()
new = open('new-1.dat','w+')
newScores = calculateMinMax(best[0],aaFreqDict2, freqDict2, mapDict, windowSize)
for k in range(0,len(newScores)):
	new.write("\n"+str(k)+" "+str(newScores[k]))
new.close()
