function [S, AllGoodNeuronIdentity]=getGoodNeurons(S,Trials)

gcs                     =load ('GoodCells_S');
gck                     =load ('GoodCells_K');
indSession              =str2num(Trials(1).blockName(3:5));

GoodCells               =[gcs.GoodCells_S; gck.GoodCells_K];
AllGoodNeuronIdentity   =GoodCells(GoodCells(:,1)==indSession,2:3);
S.AllGoodNeuronIdentity =AllGoodNeuronIdentity;

