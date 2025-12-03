function [MA,MB,MC]=date_generation(NT) 
MAdata=struct();
MBdata=struct();
MCdata=struct();

xT=linspace(0,2*pi,NT);a=(-cos(xT)+1.5);b=(cos(xT)+1.5);
MAdata.PV=1*a;
MBdata.PV=1*a;
MCdata.PV=1*a;


MAdata.D_it=10*(cos(xT)+1.5);
MBdata.D_it=10*(cos(xT+pi/6)+1.5);
MCdata.D_it=10*(cos(xT+pi/3)+1.5);

MAdata.Pmax=40;
MBdata.Pmax=40;
MCdata.Pmax=40;

MAdata.battery_capacity=40;
MBdata.battery_capacity=40;
MCdata.battery_capacity=40;

MA = Micro_Community(MAdata);
MB = Micro_Community(MBdata);
MC = Micro_Community(MCdata);
end