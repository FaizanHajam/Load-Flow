function [bus type V del P Q Qmin Qmax]=busdatafn(busdata,BMVA)
%Input the base MVA
bus=busdata(:,1);
type=busdata(:,2);
V=busdata(:,3);
del=busdata(:,4);
Pg=busdata(:,5)/BMVA;
Qg=busdata(:,6)/BMVA;
Pl=busdata(:,7)/BMVA;
Ql=busdata(:,8)/BMVA;
Qmin=busdata(:,9)/BMVA;
Qmax=busdata(:,10)/BMVA;
P=Pg-Pl;
Q=Qg-Ql;