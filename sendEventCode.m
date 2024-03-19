function sendEventCode(object,port,event)
 io64(object,port,event); %% sends event code
 WaitSecs(0.004); %% wait a few ms 
 io64(object,port,0); %% set port to zero
end