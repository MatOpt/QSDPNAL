
   function [fname,fd] = problems_sqsdp   
%%
%% theta: random
%% idx1 = [1:14];   
%%
   fname{1,1} = 'theta4';      
   fname{2,1} = 'theta42';    
   fname{3,1} = 'theta6';     
   fname{4,1} = 'theta62';    
   fname{5,1} = 'theta8';     
   fname{6,1} = 'theta82';    
   fname{7,1} = 'theta83';    
   fname{8,1} = 'theta10';    
   fname{9,1} = 'theta102';   
   fname{10,1} = 'theta103';  
   fname{11,1} = 'theta104';  
   fname{12,1} = 'theta12';   
   fname{13,1} = 'theta123';     
   fname{14,1} = 'theta162'; 
   fd(1:20) = ones(1,20);  
%%
%% theta: Dimacs 
%% idx2 = [21:38, 41:45, 51:54, 61:72];
%%
   fname{21,1} = 'MANN-a27';      
   fname{22,1} = 'johnson8-4-4';  
   fname{23,1} = 'johnson16-2-4'; 
   fname{24,1} = 'san200-0.7-1';  
   fname{25,1} = 'sanr200-0.7'; 
   fname{26,1} = 'c-fat200-1';    
   fname{27,1} = 'hamming-6-4';   
   fname{28,1} = 'hamming-8-4';   
   fname{29,1} = 'hamming-9-8';   
   fname{30,1} = 'hamming-10-2';  
   fname{31,1} = 'hamming-7-5-6'; 
   fname{32,1} = 'hamming-8-3-4'; 
   fname{33,1} = 'hamming-9-5-6'; 
   fname{34,1} = 'brock200-1';  
   fname{35,1} = 'brock200-4';  
   fname{36,1} = 'brock400-1';  
   fname{37,1} = 'keller4';     
   fname{38,1} = 'p-hat300-1';  

   fname{41,1} = 'G43';  
   fname{42,1} = 'G44';  
   fname{43,1} = 'G45';  
   fname{44,1} = 'G46';  
   fname{45,1} = 'G47';  
   fname{51,1} = 'G51'; 
   fname{52,1} = 'G52'; 
   fname{53,1} = 'G53'; 
   fname{54,1} = 'G54'; 

   fd(21:60) = 2*ones(1,40); 
%%
%% theta: Borchers
%% idx3 = [101:128];
%%
   fname{101,1} = '1dc.64';      fname{101,2} = 10;    
   fname{102,1} = '1et.64';      fname{102,2} = 18;    
   fname{103,1} = '1tc.64';      fname{103,2} = 20;    
   fname{104,1} = '1dc.128';     fname{104,2} = 16;    
   fname{105,1} = '1et.128';     fname{105,2} = 28;    
   fname{106,1} = '1tc.128';     fname{106,2} = 38;    
   fname{107,1} = '1zc.128';     fname{107,2} = 18;    
   fname{108,1} = []; '2dc.128';      fname{108,2} = 5;    
   fname{109,1} = '1dc.256';     fname{109,2} = 30;    
   fname{110,1} = '1et.256';     fname{110,2} = 50;    
   fname{111,1} = '1tc.256';     fname{111,2} = 63;    
   fname{112,1} = '1zc.256';     fname{112,2} = 36;    
   fname{113,1} = []; '2dc.256';      fname{113,2} = 7;    
   fname{114,1} = '1dc.512';     fname{114,2} = 52;    
   fname{115,1} = '1et.512';     fname{115,2} = 100;    
   fname{116,1} = '1tc.512';     fname{116,2} = 110;    
   fname{117,1} = '2dc.512';     fname{117,2} = 11;    
   fname{118,1} = '1zc.512';     fname{118,2} = 62;    
   fname{119,1} = '1dc.1024';    fname{119,2} = 94;    
   fname{120,1} = '1et.1024';    fname{120,2} = 171;    
   fname{121,1} = '1tc.1024';    fname{121,2} = 196;    
   fname{122,1} = '1zc.1024';    fname{122,2} = 112; fname{122,3} = -1;     
   fname{123,1} = '2dc.1024';    fname{123,2} = 16; fname{123,3} = -1;    
   fname{124,1} = '1dc.2048';
   fname{125,1} = '1et.2048';
   fname{126,1} = '1tc.2048';
   fname{127,1} = '1zc.2048';
   fname{128,1} = '2dc.2048';
   fname{129,1} = []; '1zc.4096';
   fd(101:130) = 3*ones(1,30); 
%%
%% FAP: 
%%
   fname{161,1} = 'fap01'; 
   fname{162,1} = 'fap02'; 
   fname{163,1} = 'fap03'; 
   fname{164,1} = 'fap04'; 
   fname{165,1} = 'fap05'; 
   fname{166,1} = 'fap06'; 
   fname{167,1} = 'fap07'; 
   fname{168,1} = 'fap08'; 
   fname{169,1} = 'fap09'; 
   fname{170,1} = 'fap10'; 
   fname{171,1} = 'fap11'; 
   fname{172,1} = 'fap12'; 
   fname{173,1} = 'fap25'; 
   fname{174,1} = 'fap36'; 
   fd(161:180) = 4*ones(1,20); 
%%
%% QAP
%%
   fname{201,1} = 'bur26a'; fname{201,2} = 5426670; 
   fname{202,1} = 'bur26b'; fname{202,2} = 3817852; 
   fname{203,1} = 'bur26c'; fname{203,2} = 5426795; 
   fname{204,1} = 'bur26d'; fname{204,2} = 3821225; 
   fname{205,1} = 'bur26e'; fname{205,2} = 5386879;
   fname{206,1} = 'bur26f'; fname{206,2} = 3782044;
   fname{207,1} = 'bur26g'; fname{207,2} = 10117172;
   fname{208,1} = 'bur26h'; fname{208,2} = 7098658;
   fname{209,1} = 'chr12a'; fname{209,2} = 9552;
   fname{210,1} = 'chr12b'; fname{210,2} = 9742;
   fname{211,1} = 'chr12c'; fname{211,2} = 11156;
   fname{212,1} = 'chr15a'; fname{212,2} = 9896;
   fname{213,1} = 'chr15b'; fname{213,2} = 7990;
   fname{214,1} = 'chr15c'; fname{214,2} = 9504;
   fname{215,1} = 'chr18a'; fname{215,2} = 11098;
   fname{216,1} = 'chr18b'; fname{216,2} = 1534;
   fname{217,1} = 'chr20a'; fname{217,2} = 2192;
   fname{218,1} = 'chr20b'; fname{218,2} = 2298;
   fname{219,1} = 'chr20c'; fname{219,2} = 14142;
   fname{220,1} = 'chr22a'; fname{220,2} = 6156;
   fname{221,1} = 'chr22b'; fname{221,2} = 6194;
   fname{222,1} = 'chr25a'; fname{222,2} = 3796;
   fname{223,1} = 'els19';  fname{223,2} = 17212548; 
   fname{224,1} = 'esc16a'; fname{224,2} = 68;
   fname{225,1} = 'esc16b'; fname{225,2} = 292;
   fname{226,1} = 'esc16c'; fname{226,2} = 160;
   fname{227,1} = 'esc16d'; fname{227,2} = 16;
   fname{228,1} = 'esc16e'; fname{228,2} = 28;
   fname{229,1} = []; 
   fname{230,1} = 'esc16g'; fname{230,2} = 26;
   fname{231,1} = 'esc16h'; fname{231,2} = 996;
   fname{232,1} = 'esc16i'; fname{232,2} = 14;
   fname{233,1} = 'esc16j'; fname{233,2} = 8;
   fname{234,1} = 'esc32a'; fname{234,2} = 130; fname{234,3} = -1; 
   fname{235,1} = 'esc32b'; fname{235,2} = 168; fname{235,3} = -1; 
   fname{236,1} = 'esc32c'; fname{236,2} = 642; fname{236,3} = -1; 
   fname{237,1} = 'esc32d'; fname{237,2} = 200; fname{237,3} = -1; 
   fname{238,1} = 'esc32e'; fname{238,2} = 2;
   fname{239,1} = 'esc32f'; fname{239,2} = 2;
   fname{240,1} = 'esc32g'; fname{240,2} = 6;
   fname{241,1} = 'esc32h'; fname{241,2} = 438; fname{241,3} = -1; 
   fname{242,1} = 'had12';  fname{242,2} = 1652;
   fname{243,1} = 'had14';  fname{243,2} = 2724;
   fname{244,1} = 'had16';  fname{244,2} = 3720;
   fname{245,1} = 'had18';  fname{245,2} = 5358;
   fname{246,1} = 'had20';  fname{246,2} = 6922;
   fname{247,1} = 'kra30a'; fname{247,2} = 88900; 
   fname{248,1} = 'kra30b'; fname{248,2} = 91420;
   fname{249,1} = 'kra32';  fname{249,2} = 88900;
   fname{250,1} = 'lipa20a'; fname{250,2} = 3683;
   fname{251,1} = 'lipa20b'; fname{251,2} = 27076;
   fname{252,1} = 'lipa30a'; fname{252,2} = 13178;
   fname{253,1} = 'lipa30b'; fname{253,2} = 151426;
   fname{254,1} = 'lipa40a'; fname{254,2} = 31538;
   fname{255,1} = 'lipa40b'; fname{255,2} = 476581;
   fname{256,1} = 'nug12';  fname{256,2} = 578;
   fname{257,1} = 'nug14';  fname{257,2} = 1014;
   fname{258,1} = 'nug15';  fname{258,2} = 1150;
   fname{259,1} = 'nug16a'; fname{259,2} = 1610;
   fname{260,1} = 'nug16b'; fname{260,2} = 1240;
   fname{261,1} = 'nug17';  fname{261,2} = 1732;
   fname{262,1} = 'nug18';  fname{262,2} = 1930;
   fname{263,1} = 'nug20';  fname{263,2} = 2570;
   fname{264,1} = 'nug21';  fname{264,2} = 2438;
   fname{265,1} = 'nug22';  fname{265,2} = 3596;
   fname{266,1} = 'nug24';  fname{266,2} = 3488;
   fname{267,1} = 'nug25';  fname{267,2} = 3744;
   fname{268,1} = 'nug27';  fname{268,2} = 5234;
   fname{269,1} = 'nug28';  fname{269,2} = 5166;
   fname{270,1} = 'nug30';  fname{270,2} = 6124;
   fname{271,1} = 'rou12';  fname{271,2} = 235528;
   fname{272,1} = 'rou15';  fname{272,2} = 354210;
   fname{273,1} = 'rou20';  fname{273,2} = 725522;
   fname{274,1} = 'scr12';  fname{274,2} = 31410;
   fname{275,1} = 'scr15';  fname{275,2} = 51140;
   fname{276,1} = 'scr20';  fname{276,2} = 110030;
   fname{277,1} = 'ste36a'; fname{277,2} = 9526;
   fname{278,1} = 'ste36b'; fname{278,2} = 15852;
   fname{279,1} = 'ste36c'; fname{279,2} = 8239110; 
   fname{280,1} = 'tai12a'; fname{280,2} = 224416; 
   fname{281,1} = 'tai12b'; fname{281,2} = 39464925;
   fname{282,1} = 'tai15a'; fname{282,2} = 388214;
   fname{283,1} = 'tai15b'; fname{283,2} = 51765268;
   fname{284,1} = 'tai17a'; fname{284,2} = 491812;
   fname{285,1} = 'tai20a'; fname{285,2} = 703482;
   fname{286,1} = 'tai20b'; fname{286,2} = 122455319;
   fname{287,1} = 'tai25a'; fname{287,2} = 1167256;
   fname{288,1} = 'tai25b'; fname{288,2} = 344355646;
   fname{289,1} = 'tai30a'; fname{289,2} = 1818146; fname{289,3} = -1;  
   fname{290,1} = 'tai30b'; fname{290,2} = 637117113; 
   fname{291,1} = 'tai35a'; fname{291,2} = 2422002;
   fname{292,1} = 'tai35b'; fname{292,2} = 283315445;
   fname{293,1} = 'tai40a'; fname{293,2} = 3139370;
   fname{294,1} = 'tai40b'; fname{294,2} = 637250948;
   fname{295,1} = 'tho30';  fname{295,2} = 149936;
   fname{296,1} = 'tho40';  fname{296,2} = 240516; fname{296,3} = -1;
   fname{297,1} = [];    
%%
   fname{298,1} = 'tai10a';[];  
   fname{299,1} = 'tai10b';[];  
   fname{300,1} = [];'esc128'; 
   fname{301,1} = []; 'esc64a';
   fname{302,1} = 'lipa50a'; 
   fname{303,1} = 'lipa50b'; 
   fname{304,1} = 'lipa60a'; 
   fname{305,1} = 'lipa60b'; 
   fname{306,1} = 'lipa70a';
   fname{307,1} = 'lipa70b'; 
   fname{308,1} = 'lipa80a'; 
   fname{309,1} = 'lipa80b'; 
   fname{310,1} = 'lipa90a'; 
   fname{311,1} = 'lipa90b';
   fname{312,1} = 'sko100a'; 
   fname{313,1} = 'sko100b'; 
   fname{314,1} = 'sko100c'; 
   fname{315,1} = 'sko100d'; 
   fname{316,1} = 'sko100e';
   fname{317,1} = 'sko100f'; 
   fname{318,1} = 'sko42'; 
   fname{319,1} = 'sko49'; 
   fname{320,1} = 'sko56'; 
   fname{321,1} = 'sko64'; 
   fname{322,1} = 'sko72'; 
   fname{323,1} = 'sko81'; 
   fname{324,1} = 'sko90';  
   fname{325,1} = 'tai100a'; 
   fname{326,1} = 'tai100b'; 
   fname{327,1} = 'tai50a'; 
   fname{328,1} = 'tai50b'; 
   fname{329,1} = 'tai60a'; 
   fname{330,1} = 'tai60b'; 
   fname{331,1} = 'tai64c'; 
   fname{332,1} = 'tai80a'; 
   fname{333,1} = 'tai80b';
   fname{334,1} = 'tai256c'; 
   fname{335,1} = 'tai150b'; 
   fname{336,1} = 'tho150';
   fname{337,1} = 'wil100'; 
   fname{338,1} = 'wil50';

   fd(201:340) = 5*ones(1,140); 
%%
%% NCM
%%
   fname{401,1} = 'NCM1n200H1'; 
   fname{402,1} = 'NCM1n200H2'; 
   fname{403,1} = 'NCM1n400H1'; 
   fname{404,1} = 'NCM1n400H2';
   fname{405,1} = 'NCM1n800H1';
   fname{406,1} = 'NCM1n800H2';
  
   fd(401:410) = 6*ones(1,10); 

%%
%% sparse maximum eigenvalue
%%
   fname{411,1} = 'spmaxeig4n200H1'; 
   fname{412,1} = 'spmaxeig5n200H2'; 
   fname{413,1} = 'spmaxeig1n400H1'; 
   fname{414,1} = 'spmaxeig4n400H2';
   fname{415,1} = 'spmaxeig1n800H1';
   fname{416,1} = 'spmaxeig4n800H2';
  
   fd(411:420) = 6.1*ones(1,10); 
%%
%% molecular conformation
%%
   fname{430,1} = 'PTQ30n2';

   fd(431:440)  = 6.5*ones(1,10); 
%%
%% fukuda: quantum chemistry
%%
   fname{451,1} = 'BeH_2Sigma+_STO-6GN5r12g1T2';
   fname{452,1} = 'BH_1Sigma+_STO-6GN6r12g1T2';
   fname{453,1} = 'BH2_2A1_STO-6GN7r14g1T2';
   fname{454,1} = 'BH+_2Sigma+_STO-6GN5r12g1T2';
   fname{455,1} = 'CH+_1Sigma+_STO-6GN6r12g1T2';
   fname{456,1} = 'CH2_1A1_STO-6GN8r14g1T2';
   fname{457,1} = 'CH2_3B1_STO-6GN8r14g1T2';
   fname{458,1} = 'CH_2Pi_STO-6GN7r12g1T2';
   fname{459,1} = 'CH-_3Sigma-_STO-6GN8r12g1T2';
   fname{460,1} = 'H2O_1A1_STO-6GN10r14g1T2';
   fname{461,1} = 'H2O+_2B1_STO-6GN9r14g1T2';
   fname{462,1} = 'HF_1Sigma+_STO-6GN10r12g1T2_5';
   fname{463,1} = 'HF+_2Pi_STO-6GN9r12g1T2';
   fname{464,1} = 'LiH_1Sigma+_STO-6GN4r12g1T2';
   fname{465,1} = 'NH2_2B1_STO-6GN9r14g1T2';
   fname{466,1} = 'NH+_2Pi_STO-6GN7r12g1T2';
   fname{467,1} = 'NH-_2Pi_STO-6GN9r12g1T2';
   fname{468,1} = 'NH_3Sigma-_STO-6GN8r12g1T2';
   fname{469,1} = 'OH-_1Sigma+_STO-6GN10r12g1T2_5';
   fname{470,1} = 'OH_2Pi_STO-6GN9r12g1T2';
   fname{471,1} = 'OH+_3Sigma-_STO-6GN8r12g1T2';
   fname{472,1} = 'H3O+_1-A1_STO-6GN10r16g1T2_5';
   fname{473,1} = 'NH3_1-A1_STO-6GN10r16g1T2_5';

   fd([451:480]) = 6.7*ones(1,30);    
%%
%% sparse random SDPs of Rendl
%%
   fname{501,1} = 'Rn3m20p3'; 
   fname{502,1} = 'Rn3m25p3'; 
   fname{503,1} = 'Rn3m10p4'; 

   fname{504,1} = 'Rn4m30p3'; 
   fname{505,1} = 'Rn4m40p3'; 
   fname{506,1} = 'Rn4m15p4'; 

   fname{507,1} = 'Rn5m30p3'; 
   fname{508,1} = 'Rn5m40p3'; 
   fname{509,1} = 'Rn5m50p3'; 
   fname{510,1} = 'Rn5m20p4'; 

   fname{511,1} = 'Rn6m40p3'; 
   fname{512,1} = 'Rn6m50p3'; 
   fname{513,1} = 'Rn6m60p3'; 
   fname{514,1} = 'Rn6m20p4'; 

   fname{515,1} = 'Rn7m50p3'; 
   fname{516,1} = 'Rn7m70p3'; 
   fname{517,1} = 'Rn8m70p3'; 
   fname{518,1} = 'Rn8m100p3'; 

   fname{530,1} = []; 'Rn2m10p3'; 

   fd(501:530) = 7*ones(1,30); 

%%
%% binary quadratic programming
%%
    fname{601,1} = 'be100.1';   fname{601,2} = -19412;
    fname{602,1} = 'be100.2';  fname{602,2} = -17290;
    fname{603,1} = 'be100.3';  fname{603,2} = -17565;   
    fname{604,1} = 'be100.4';  fname{604,2} = -19125;   
    fname{605,1} = 'be100.5';  fname{605,2} = -15868;
    fname{606,1} = 'be100.6';  fname{606,2} = -17368;   
    fname{607,1} = 'be100.7';  fname{607,2} = -18629;   
    fname{608,1} = 'be100.8';  fname{608,2} = -18649;
    fname{609,1} = 'be100.9';  fname{609,2} = -13294;   
    fname{610,1} = 'be100.10'; fname{610,2} = -15352;
    fname{611,1} = 'be120.3.1'; fname{611,2} = -13067;
    fname{612,1} = 'be120.3.2'; fname{612,2} = -13046;  
    fname{613,1} = 'be120.3.3'; fname{613,2} = -12418;    
    fname{614,1} = 'be120.3.4'; fname{614,2} = -13867;  
    fname{615,1} = 'be120.3.5'; fname{615,2} = -11403;    
    fname{616,1} = 'be120.3.6'; fname{616,2} = -12915;    
    fname{617,1} = 'be120.3.7'; fname{617,2} = -14068;  
    fname{618,1} = 'be120.3.8'; fname{618,2} = -14701;    
    fname{619,1} = 'be120.3.9'; fname{619,2} = -10458;  
    fname{620,1} = 'be120.3.10'; fname{620,2} = -12201;      
    fname{621,1} = 'be120.8.1';  fname{621,2} = -18691;    
    fname{622,1} = 'be120.8.2';  fname{622,2} = -18827;     
    fname{623,1} = 'be120.8.3';  fname{623,2} = -19302;    
    fname{624,1} = 'be120.8.4';  fname{624,2} = -20765;     
    fname{625,1} = 'be120.8.5';  fname{625,2} = -20417;     
    fname{626,1} = 'be120.8.6';  fname{626,2} = -18482;    
    fname{627,1} = 'be120.8.7';  fname{627,2} = -22194;     
    fname{628,1} = 'be120.8.8';  fname{628,2} = -19534;     
    fname{629,1} = 'be120.8.9';  fname{629,2} = -18195;    
    fname{630,1} = 'be120.8.10'; fname{630,2} = -19049;    
    fname{631,1} = 'be150.3.1';  fname{631,2} = -18889;    
    fname{632,1} = 'be150.3.2';  fname{632,2} = -17816;    
    fname{633,1} = 'be150.3.3';  fname{633,2} = -17314;    
    fname{634,1} = 'be150.3.4';  fname{634,2} = -19884;     
    fname{635,1} = 'be150.3.5';  fname{635,2} = -16817;    
    fname{636,1} = 'be150.3.6';  fname{636,2} = -16780;      
    fname{637,1} = 'be150.3.7';  fname{637,2} = -18001;    
    fname{638,1} = 'be150.3.8';  fname{638,2} = -18303;    
    fname{639,1} = 'be150.3.9';  fname{639,2} = -12838;    
    fname{640,1} = 'be150.3.10'; fname{640,2} = -17963;        
    fname{641,1} = 'be150.8.1';  fname{641,2} = -27089;        
    fname{642,1} = 'be150.8.2';  fname{642,2} = -26779;         
    fname{643,1} = 'be150.8.3';  fname{643,2} = -29438;         
    fname{644,1} = 'be150.8.4';  fname{644,2} = -26911;        
    fname{645,1} = 'be150.8.5';  fname{645,2} = -28017;         
    fname{646,1} = 'be150.8.6';  fname{646,2} = -29221;         
    fname{647,1} = 'be150.8.7';  fname{647,2} = -31209;        
    fname{648,1} = 'be150.8.8';  fname{648,2} = -29730;         
    fname{649,1} = 'be150.8.9';  fname{649,2} = -25388;         
    fname{650,1} = 'be150.8.10'; fname{650,2} = -28374;         
    fname{651,1} = 'be200.3.1';  fname{651,2} = -25453;   
    fname{652,1} = 'be200.3.2';  fname{652,2} = -25027;   
    fname{653,1} = 'be200.3.3';  fname{653,2} = -28023;   
    fname{654,1} = 'be200.3.4';  fname{654,2} = -27434;   
    fname{655,1} = 'be200.3.5';  fname{655,2} = -26355;   
    fname{656,1} = 'be200.3.6';  fname{656,2} = -26146;   
    fname{657,1} = 'be200.3.7';  fname{657,2} = -30483;   
    fname{658,1} = 'be200.3.8';  fname{658,2} = -27355;   
    fname{659,1} = 'be200.3.9';  fname{659,2} = -24683;   
    fname{660,1} = 'be200.3.10'; fname{660,2} = -23842;   
    fname{661,1} = 'be200.8.1';  fname{661,2} = -48534;    
    fname{662,1} = 'be200.8.2';  fname{662,2} = -40821;   
    fname{663,1} = 'be200.8.3';  fname{663,2} = -43207;   
    fname{664,1} = 'be200.8.4';  fname{664,2} = -43757;   
    fname{665,1} = 'be200.8.5';  fname{665,2} = -41482;   
    fname{666,1} = 'be200.8.6';  fname{666,2} = -49492;   
    fname{667,1} = 'be200.8.7';  fname{667,2} = -46828;   
    fname{668,1} = 'be200.8.8';  fname{668,2} = -44502;   
    fname{669,1} = 'be200.8.9';  fname{669,2} = -43241;   
    fname{670,1} = 'be200.8.10'; fname{670,2} = -42832;   
    fname{671,1} = 'be250.1';    fname{671,2} = -24076;   
    fname{672,1} = 'be250.2';    fname{672,2} = -22540;   
    fname{673,1} = 'be250.3';    fname{673,2} = -22923;   
    fname{674,1} = 'be250.4';    fname{674,2} = -24649;   
    fname{675,1} = 'be250.5';    fname{675,2} = -21057;   
    fname{676,1} = 'be250.6';    fname{676,2} = -22735;   
    fname{677,1} = 'be250.7';    fname{677,2} = -24095;   
    fname{678,1} = 'be250.8';    fname{678,2} = -23801;   
    fname{679,1} = 'be250.9';    fname{679,2} = -20051;   
    fname{680,1} = 'be250.10';   fname{680,2} = -23159;   
    fname{681,1} = 'bqp50-1';    fname{681,2} = -2098;        
    fname{682,1} = 'bqp50-2';    fname{682,2} = -3702;                
    fname{683,1} = 'bqp50-3';    fname{683,2} = -4626;                  
    fname{684,1} = 'bqp50-4';    fname{684,2} = -3544;                
    fname{685,1} = 'bqp50-5';    fname{685,2} = -4012;                  
    fname{686,1} = 'bqp50-6';    fname{686,2} = -3693;                 
    fname{687,1} = 'bqp50-7';    fname{687,2} = -4520;             
    fname{688,1} = 'bqp50-8';    fname{688,2} = -4216;                 
    fname{689,1} = 'bqp50-9';    fname{689,2} = -3780;        
    fname{690,1} = 'bqp50-10';   fname{690,2} = -3507;                      
    fname{691,1} = 'bqp100-1';   fname{691,2} = -7970;                        
    fname{692,1} = 'bqp100-2';   fname{692,2} = -11036;       
    fname{693,1} = 'bqp100-3';   fname{693,2} = -12723;            
    fname{694,1} = 'bqp100-4';   fname{694,2} = -10368;          
    fname{695,1} = 'bqp100-5';   fname{695,2} = -9083;           
    fname{696,1} = 'bqp100-6';   fname{696,2} = -10210;          
    fname{697,1} = 'bqp100-7';   fname{697,2} = -10125;       
    fname{698,1} = 'bqp100-8';   fname{698,2} = -11435;          
    fname{699,1} = 'bqp100-9';   fname{699,2} = -11455;           
    fname{700,1} = 'bqp100-10';  fname{700,2} = -12565;           
    fname{701,1} = 'bqp250-1';   fname{701,2} = -45607;        
    fname{702,1} = 'bqp250-2';   fname{702,2} = -44810;             
    fname{703,1} = 'bqp250-3';   fname{703,2} = -49037;         
    fname{704,1} = 'bqp250-4';   fname{704,2} = -41274;            
    fname{705,1} = 'bqp250-5';   fname{705,2} = -47961;             
    fname{706,1} = 'bqp250-6';   fname{706,2} = -41014;            
    fname{707,1} = 'bqp250-7';   fname{707,2} = -46757;             
    fname{708,1} = 'bqp250-8';   fname{708,2} = -35726;           
    fname{709,1} = 'bqp250-9';   fname{709,2} = -48916;        
    fname{710,1} = 'bqp250-10';  fname{710,2} = -40442;           
    fname{711,1} = 'bqp500-1';   fname{711,2} = -116586;                 
    fname{712,1} = 'bqp500-2';   fname{712,2} = -128223;                      
    fname{713,1} = 'bqp500-3';   fname{713,2} = -130812;                       
    fname{714,1} = 'bqp500-4';   fname{714,2} = -130097;               
    fname{715,1} = 'bqp500-5';   fname{715,2} = -125487;                       
    fname{716,1} = 'bqp500-6';   fname{716,2} = -121772;           
    fname{717,1} = 'bqp500-7';   fname{717,2} = -122201;  
    fname{718,1} = 'bqp500-8';   fname{718,2} = -123559;          
    fname{719,1} = 'bqp500-9';   fname{719,2} = -120798;           
    fname{720,1} = 'bqp500-10';  fname{720,2} = -130619;      
    fname{721,1} = 'gka1a';      fname{721,2} = -3414;      
    fname{722,1} = 'gka2a';      fname{722,2} = -6063;      
    fname{723,1} = 'gka3a';      fname{723,2} = -6037;      
    fname{724,1} = 'gka4a';      fname{724,2} = -8598;       
    fname{725,1} = 'gka5a';      fname{725,2} = -5737;      
    fname{726,1} = 'gka6a';      fname{726,2} = -3980;      
    fname{727,1} = 'gka7a';      fname{727,2} = -4541;      
    fname{728,1} = 'gka8a';      fname{728,2} = -11109;       
    fname{729,1} = 'gka1b';  fname{729,2} = -133;         
    fname{730,1} = 'gka2b';  fname{730,2} = -121;         
    fname{731,1} = 'gka3b';  fname{731,2} = -118;              
    fname{732,1} = 'gka4b';  fname{732,2} = -129;             
    fname{733,1} = 'gka5b';  fname{733,2} = -150;             
    fname{734,1} = 'gka6b';  fname{734,2} = -146;             
    fname{735,1} = 'gka7b';  fname{735,2} = -160;             
    fname{736,1} = 'gka8b';  fname{736,2} = -145;             
    fname{737,1} = 'gka9b';  fname{737,2} = -137;             
    fname{738,1} = 'gka10b'; fname{738,2} = -154;              
    fname{739,1} = 'gka1c';  fname{739,2} = -5058;         
    fname{740,1} = 'gka2c';  fname{740,2} = -6213;               
    fname{741,1} = 'gka3c';  fname{741,2} = -6665;             
    fname{742,1} = 'gka4c';  fname{742,2} = -7398;                 
    fname{743,1} = 'gka5c';  fname{743,2} = -7362;                 
    fname{744,1} = 'gka6c';  fname{744,2} = -5824;             
    fname{745,1} = 'gka7c';  fname{745,2} = -7225;                 
    fname{746,1} = 'gka1d';  fname{746,2} = -6333;                  
    fname{747,1} = 'gka2d';  fname{747,2} = -6579;                      
    fname{748,1} = 'gka3d';  fname{748,2} = -9261;                      
    fname{749,1} = 'gka4d';  fname{749,2} = -10727;                      
    fname{750,1} = 'gka5d';  fname{750,2} = -11626;                  
    fname{751,1} = 'gka6d';  fname{751,2} = -14207;                       
    fname{752,1} = 'gka7d';  fname{752,2} = -14476;                  
    fname{753,1} = 'gka8d';  fname{753,2} = -16352;                      
    fname{754,1} = 'gka9d';  fname{754,2} = -15656;                  
    fname{755,1} = 'gka10d'; fname{755,2} = -19102;                      
    fname{756,1} = 'gka1e';  fname{756,2} = -16464;                      
    fname{757,1} = 'gka2e';  fname{757,2} = -23395;                      
    fname{758,1} = 'gka3e';  fname{758,2} = -25243;                      
    fname{759,1} = 'gka4e';  fname{759,2} = -35594;                  
    fname{760,1} = 'gka5e';  fname{760,2} = -35154;                       
    fname{761,1} = 'gka1f';  fname{761,2} = -61194;   %% non-optimal 
    fname{762,1} = 'gka2f';  fname{762,2} = -100161;  %% non-optimal 
    fname{763,1} = 'gka3f';  fname{763,2} = -138035;  %% non-optimal
    fname{764,1} = 'gka4f';  fname{764,2} = -172771;  %% non-optimal  
    fname{765,1} = 'gka5f';  fname{765,2} = -190507;  %% non-optimal 

    fd(601:780) = 8*ones(1,180); 
%%
%% random symmetric matrix completion
%%
   fname{801,1} = 'Sn0500r010d05e0';
   fname{802,1} = 'Sn0500r010d05e5';
   fname{803,1} = 'Sn1000r020d05e0';
   fname{804,1} = 'Sn1000r020d05e5';
   fname{805,1} = 'Sn2000r040d05e0';
   fname{806,1} = 'Sn2000r040d05e5';
   fname{807,1} = 'Sn1000r050d05e0';
   fname{808,1} = 'Sn1000r050d05e5';
   fname{809,1} = 'Sn2000r100d05e0';
   fname{810,1} = 'Sn2000r100d05e5';
   fd(801:810) = 9*ones(1,10); 
%%
%% Molecular conformation
%%

   fname{901,1} = '1GM2';
   fname{902,1} = '1PTQ';
   fname{903,1} = '1AX8';
   fname{904,1} = '1F39';
   fname{905,1} = '1TIMa';
   fname{906,1} = '1RGS';
   fname{907,1} = '1KDH';

   fd(901:910) = 10*ones(1,10); 
%%
%% maxcut
%%
%    fname{1001,1} = 'G43';  
%    fname{1002,1} = 'G44';  
%    fname{1003,1} = 'G45';  
%    fname{1004,1} = 'G46';  
%    fname{1005,1} = 'G47';  
%    fname{1006,1} = 'G51'; 
%    fname{1007,1} = 'G52'; 
%    fname{1008,1} = 'G53'; 
%    fname{1009,1} = 'G54'; 
%    fname{1010,1} = 'G32'; 

   fd(1001:1020) = 11*ones(1,20); 
%%
%% Shirai
%%
   fname{1021,1} = 'adet3_2.0_N=2_';
   fname{1022,1} = 'adet3_2.0_N=3_';
   fname{1023,1} = 'adet3_2.0_N=4_';
   fname{1024,1} = 'adet4_2.0_N=3_';
   fname{1025,1} = 'adet5_2.0_N=4';
   fname{1026,1} = 'adet6_2.0_N=3';
   
   fd(1021:1030) = 12*ones(1,10);
%%**************************************************************

%% relaxed clustering problem
%%
   
   fname{1041,1} = 'soybean-small-2';
   fname{1042,1} = 'soybean-small-3';
   fname{1043,1} = 'soybean-small-4';
   fname{1044,1} = 'soybean-small-5';
   fname{1045,1} = 'soybean-small-6';
   fname{1046,1} = 'soybean-small-7';
   fname{1047,1} = 'soybean-small-8';
   fname{1048,1} = 'soybean-small-9';
   fname{1049,1} = 'soybean-small-10';
   fname{1050,1} = 'soybean-small-11';
   
   fname{1051,1} = 'soybean-large-2';
   fname{1052,1} = 'soybean-large-3';
   fname{1053,1} = 'soybean-large-4';
   fname{1054,1} = 'soybean-large-5';
   fname{1055,1} = 'soybean-large-6';
   fname{1056,1} = 'soybean-large-7';
   fname{1057,1} = 'soybean-large-8';
   fname{1058,1} = 'soybean-large-9';
   fname{1059,1} = 'soybean-large-10';
   fname{1060,1} = 'soybean-large-11';
   
   fname{1061,1} = 'spambase-small-2';
   fname{1062,1} = 'spambase-small-3';
   fname{1063,1} = 'spambase-small-4';
   fname{1064,1} = 'spambase-small-5';
   fname{1065,1} = 'spambase-small-6';
   fname{1066,1} = 'spambase-small-7';
   fname{1067,1} = 'spambase-small-8';
   fname{1068,1} = 'spambase-small-9';
   fname{1069,1} = 'spambase-small-10';
   fname{1070,1} = 'spambase-small-11';

   fname{1071,1} = 'spambase-medium-2';
   fname{1072,1} = 'spambase-medium-3';
   fname{1073,1} = 'spambase-medium-4';
   fname{1074,1} = 'spambase-medium-5';
   fname{1075,1} = 'spambase-medium-6';
   fname{1076,1} = 'spambase-medium-7';
   fname{1077,1} = 'spambase-medium-8';
   fname{1078,1} = 'spambase-medium-9';
   fname{1079,1} = 'spambase-medium-10';
   fname{1080,1} = 'spambase-medium-11';
   
   fname{1081,1} = 'spambase-large-2';
   fname{1082,1} = 'spambase-large-3';
   fname{1083,1} = 'spambase-large-4';
   fname{1084,1} = 'spambase-large-5';
   fname{1085,1} = 'spambase-large-6';
   fname{1086,1} = 'spambase-large-7';
   fname{1087,1} = 'spambase-large-8';
   fname{1088,1} = 'spambase-large-9';
   fname{1089,1} = 'spambase-large-10';
   fname{1090,1} = 'spambase-large-11';
   
   fname{1091,1} = 'abalone-small-2';
   fname{1092,1} = 'abalone-small-3';
   fname{1093,1} = 'abalone-small-4';
   fname{1094,1} = 'abalone-small-5';
   fname{1095,1} = 'abalone-small-6';
   fname{1096,1} = 'abalone-small-7';
   fname{1097,1} = 'abalone-small-8';
   fname{1098,1} = 'abalone-small-9';
   fname{1099,1} = 'abalone-small-10';
   fname{1100,1} = 'abalone-small-11';
   
   fname{1101,1} = 'abalone-medium-2';
   fname{1102,1} = 'abalone-medium-3';
   fname{1103,1} = 'abalone-medium-4';
   fname{1104,1} = 'abalone-medium-5';
   fname{1105,1} = 'abalone-medium-6';
   fname{1106,1} = 'abalone-medium-7';
   fname{1107,1} = 'abalone-medium-8';
   fname{1108,1} = 'abalone-medium-9';
   fname{1109,1} = 'abalone-medium-10';
   fname{1110,1} = 'abalone-medium-11';
   
   fname{1111,1} = 'abalone-large-2';
   fname{1112,1} = 'abalone-large-3';
   fname{1113,1} = 'abalone-large-4';
   fname{1114,1} = 'abalone-large-5';
   fname{1115,1} = 'abalone-large-6';
   fname{1116,1} = 'abalone-large-7';
   fname{1117,1} = 'abalone-large-8';
   fname{1118,1} = 'abalone-large-9';
   fname{1119,1} = 'abalone-large-10';
   fname{1120,1} = 'abalone-large-11';
   
   fname{1121,1} = 'segment-small-2';
   fname{1122,1} = 'segment-small-3';
   fname{1123,1} = 'segment-small-4';
   fname{1124,1} = 'segment-small-5';
   fname{1125,1} = 'segment-small-6';
   fname{1126,1} = 'segment-small-7';
   fname{1127,1} = 'segment-small-8';
   fname{1128,1} = 'segment-small-9';
   fname{1129,1} = 'segment-small-10';
   fname{1130,1} = 'segment-small-11';
   
   fname{1131,1} = 'segment-medium-2';
   fname{1132,1} = 'segment-medium-3';
   fname{1133,1} = 'segment-medium-4';
   fname{1134,1} = 'segment-medium-5';
   fname{1135,1} = 'segment-medium-6';
   fname{1136,1} = 'segment-medium-7';
   fname{1137,1} = 'segment-medium-8';
   fname{1138,1} = 'segment-medium-9';
   fname{1139,1} = 'segment-medium-10';
   fname{1140,1} = 'segment-medium-11';
   
   fname{1141,1} = 'segment-large-2';
   fname{1142,1} = 'segment-large-3';
   fname{1143,1} = 'segment-large-4';
   fname{1144,1} = 'segment-large-5';
   fname{1145,1} = 'segment-large-6';
   fname{1146,1} = 'segment-large-7';
   fname{1147,1} = 'segment-large-8';
   fname{1148,1} = 'segment-large-9';
   fname{1149,1} = 'segment-large-10';
   fname{1150,1} = 'segment-large-11';
   
   fname{1151,1} = 'housing-2'; 
   fname{1152,1} = 'housing-3';
   fname{1153,1} = 'housing-4'; 
   fname{1154,1} = 'housing-5';
   fname{1155,1} = 'housing-6'; 
   fname{1156,1} = 'housing-7';
   fname{1157,1} = 'housing-8'; 
   fname{1158,1} = 'housing-9';
   fname{1159,1} = 'housing-10'; 
   fname{1160,1} = 'housing-11';
   fd([1041:1200]) = 13*ones(1,160);
%%**************************************************************

%% Rank-1 tensor approximation
%%  
   fname{1801,1} = 'nonsym(5,4)';
   fname{1802,1} = 'nonsym(6,4)';
   fname{1803,1} = 'nonsym(7,4)';
   fname{1804,1} = 'nonsym(8,4)';
   fname{1805,1} = 'nonsym(9,4)';
   fname{1806,1} = 'nonsym(10,4)';
   fname{1807,1} = 'nonsym(11,4)';
   fname{1808,1} = 'nonsym(3,5)';
   fname{1809,1} = 'nonsym(4,5)';
   fname{1810,1} = 'nonsym(5,5)';
   fname{1811,1} = 'nonsym(6,5)';
   
   fname{1812,1} = 'sym_rd(3,20)';
   fname{1813,1} = 'sym_rd(3,25)';
   fname{1814,1} = 'sym_rd(3,30)';
   fname{1815,1} = 'sym_rd(3,35)';
   fname{1816,1} = 'sym_rd(3,40)';
   fname{1817,1} = 'sym_rd(3,45)';
   fname{1818,1} = 'sym_rd(3,50)';
   fname{1819,1} = 'sym_rd(4,20)';
   fname{1820,1} = 'sym_rd(4,25)';
   fname{1821,1} = 'sym_rd(4,30)';
   fname{1822,1} = 'sym_rd(4,35)';
   fname{1823,1} = 'sym_rd(4,40)';
   fname{1824,1} = 'sym_rd(4,45)';
   fname{1825,1} = 'sym_rd(4,50)';
   fname{1826,1} = 'sym_rd(5,5)';
   fname{1827,1} = 'sym_rd(5,10)';
   fname{1828,1} = 'sym_rd(5,15)';
   fname{1829,1} = 'sym_rd(5,20)';
   fname{1830,1} = 'sym_rd(6,5)';
   fname{1831,1} = 'sym_rd(6,10)';
   fname{1832,1} = 'sym_rd(6,15)';
   fname{1833,1} = 'sym_rd(6,20)';
   
   fname{1834,1} = 'nsym_rd([10,10,10])';
   fname{1835,1} = 'nsym_rd([15,15,15])';
   fname{1836,1} = 'nsym_rd([20,20,20])';
   fname{1837,1} = 'nsym_rd([20,25,25])';
   fname{1838,1} = 'nsym_rd([25,20,25])';
   fname{1839,1} = 'nsym_rd([25,25,20])';
   fname{1840,1} = 'nsym_rd([25,25,25])';
   fname{1841,1} = 'nsym_rd([30,30,30])';
   fname{1842,1} = 'nsym_rd([35,35,35])';
   fname{1843,1} = 'nsym_rd([40,40,40])';
   fname{1844,1} = 'nsym_rd([5,5,5,5])';
   fname{1845,1} = 'nsym_rd([6,6,6,6])';
   fname{1846,1} = 'nsym_rd([7,7,7,7])';
   fname{1847,1} = 'nsym_rd([8,8,8,8])';
   fname{1848,1} = 'nsym_rd([9,9,9,9])';
   
   
   fname{1849,1} = 'nonsym(12,4)';
   fname{1850,1} = 'nonsym(13,4)';
   fname{1851,1} = 'nonsym(14,4)';
   fname{1852,1} = 'nonsym(15,4)';
   fname{1853,1} = 'nonsym(16,4)';
   fname{1854,1} = 'nonsym(17,4)';
   fname{1855,1} = 'nonsym(7,5)';
   fname{1856,1} = 'nonsym(8,5)';
   
   fd([1801:2000]) = 14*ones(1,200);
%%**************************************************************


