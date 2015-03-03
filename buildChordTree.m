% the chord tree is built this way with priority from top to down
% --music--            --digital--     --digital difference--
% 1->3->5 maj          1,5,8            4,7
% 1->3b->5 min         1,4,8            3,7
% 1->3->5# aug         1,5,9            4,8
% 1->3b->5b dim        1,4,7            3,6
% 1->3->7 maj7         1,5,12           4,11
% 1->3b->7b min7       1,4,11           3,10
% 1->3b->7 minmaj7     1,4,12           3,11
% 1->3->7b dom7        1,5,11           4,10
% 1->3->6 maj6         1,5,10           4,9
% 1->3b->6 min6        1,4,10           3,9
% 1->4->5 sus4         1,6,8            5,7
% 1->2->5 sus2         1,3,8            2,7
% 1->3b->6b maj/3      1,4,9            3,8
% 1->4->6 maj/5        1,6,10           5,9
% 1->2->4 min/7        1,3,6            2,5
% 1->2->4# maj/7       1,3,7            2,6
% 1->1#->4 maj/7+      1,2,6            1,5
% 1->2->7b maj/2       1,3,11           2,10
% 1->3b min            1,4              3
% 1->3 maj             1,5              4
% 1->5 5               1,8              7
% 1->2 2               1,3              2
% 1->4 4               1,6              5
% 1->4# 4#             1,7              6
% 1->6b 6b             1,9              8
% 1->6 6               1,10             9
% 1->7b 7b             1,11             10
% 1->7 7               1,12             11
% default N            x,x,x            x,x
function chordtree = buildChordTree
nchordtype = 28;
chordtree = cell(2,nchordtype);

chordtree{1,1} = [3,7];
chordtree{2,1} = 'min';

chordtree{1,2} = [4,7];
chordtree{2,2} = 'maj';

chordtree{1,3} = [4,8];
chordtree{2,3} = 'aug';

chordtree{1,4} = [3,6];
chordtree{2,4} = 'dim';

chordtree{1,5} = [3,10];
chordtree{2,5} = 'min7';

chordtree{1,6} = [4,11];
chordtree{2,6} = 'maj7';

chordtree{1,7} = [3,11];
chordtree{2,7} = 'minmaj7';

chordtree{1,8} = [4,10];
chordtree{2,8} = 'dom7';

chordtree{1,9} = [3,9];
chordtree{2,9} = 'min6';

chordtree{1,10} = [4,9];
chordtree{2,10} = 'maj6';

chordtree{1,11} = [5,7];
chordtree{2,11} = 'sus4';

chordtree{1,12} = [2,7];
chordtree{2,12} = 'sus2';

chordtree{1,13} = [3,8];
chordtree{2,13} = 'maj/3';

chordtree{1,14} = [5,9];
chordtree{2,14} = 'maj/5';

chordtree{1,15} = [2,5];
chordtree{2,15} = 'min/7';

chordtree{1,16} = [2,6];
chordtree{2,16} = 'maj/7';

chordtree{1,17} = [1,5];
chordtree{2,17} = 'maj/7+';

chordtree{1,18} = [2,10];
chordtree{2,18} = 'maj/2';

chordtree{1,19} = 3;
chordtree{2,19} = 'min';

chordtree{1,20} = 4;
chordtree{2,20} = 'maj';

chordtree{1,21} = 7;
chordtree{2,21} = '5';

chordtree{1,22} = 2;
chordtree{2,22} = '2';

chordtree{1,23} = 5;
chordtree{2,23} = '4';

chordtree{1,24} = 6;
chordtree{2,24} = '4#';

chordtree{1,25} = 8;
chordtree{2,25} = '6b';

chordtree{1,26} = 9;
chordtree{2,26} = '6';

chordtree{1,27} = 10;
chordtree{2,27} = '7b';

chordtree{1,28} = 11;
chordtree{2,28} = '7';
