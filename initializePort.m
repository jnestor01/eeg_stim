% startEvent = 255;
function [object,port] =  initializePort(startEvent)

%%% initialize driver at the beginning of each experiment

object = io64;  %% initialize driver of parallel port
status = io64(object);  %% check status of driver (should be open now)
        %%% should say: 64-bit Windows
        %%% 'InpOut32a driver is open' 
        %%% [status = 0]
 port= 69632; %hex2dec('C010');% Parallel port's address in computer (check under devices/Resources Setting)
 data = io64(object,port);  %% read out data from port
 data_in=io64(object,port); %% read out in matlab agai
 
 
 %%% ------------------------------------------------- %%%
 %%% THIS SENDS CODES; ALWAYS SET PORT BACK TO 0 AFTER SENDING A CODE!!
 io64(object,port,0); 
 io64(object,port,startEvent); WaitSecs(0.004);%% sends event code
 io64(object,port,0); %% set port to zero

end
 
 