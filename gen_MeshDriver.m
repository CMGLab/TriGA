clc
clear 
close all

addpath ./mesh2D
load circle


[IEN, NODE,node,bflag] = gen_mesh(P,KV,h,option);

showMesh(node)