load fort.555;

time = fort(:,1);
x = fort(:,2);
height = fort(:,3);
extp = fort(:,4);

index = find(time==0);
figure(2);plot(x(index),height(index));
figure(3);plot(x(index),extp(index));


npts = length(index);
nptsStart = npts+1

while (nptsStart < length(time))
    nptsEnd = nptsStart + npts - 1;
    figure(2); plot(x(nptsStart:nptsEnd)/1000, height(nptsStart:nptsEnd));
    axis([-160,160,-.010,.010]);
    
    xlabel('distance (km)');
    ylabel('height (m)');
    
    figure(3); plot(x(nptsStart:nptsEnd)/1000, extp(nptsStart:nptsEnd));
   
    axis([-160,160,-1 ,5]);
    xlabel('distance (km)');
    ylabel('% overpressure'); 
    
    nptsStart = nptsEnd + 1;  % prepare for next iteration
    
    
end