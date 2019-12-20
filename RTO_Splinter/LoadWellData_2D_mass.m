function [ w1,w2 ] = LoadWellData_2D_mass(Pr1,Pr2)
%LOADWELLDATA Summary of this function goes here
%   Detailed explanation goes here

str1 = num2str(Pr1);
str2 = num2str(Pr2);

eval(['load ' 'IPR_' str1 '_2D_Splines_mass.mat;'])
eval(['load ' 'IPR_' str2 '_2D_Splines_mass.mat;'])

w1 = [];
w2 = [];
eval(['w1.FQG = QG_w' str1 ';'])
eval(['w2.FQG = QG_w' str2 ';'])

eval(['w1.FQL = QL_w' str1 ';'])
eval(['w2.FQL = QL_w' str2 ';'])

eval(['w1.FWG = WG_w' str1 ';'])
eval(['w2.FWG = WG_w' str2 ';'])

eval(['w1.FWL = WL_w' str1 ';'])
eval(['w2.FWL = WL_w' str2 ';'])

eval(['w1.FPWH = PWH_w' str1 ';'])
eval(['w2.FPWH = PWH_w' str2 ';'])

eval(['w1.FTWH = TWH_w' str1 ';'])
eval(['w2.FTWH = TWH_w' str2 ';'])


end

