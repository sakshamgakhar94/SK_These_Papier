figure(1)
plot(zvel(:,1),zvel(:,2));

peak_vel_threshold=6; % A MODULER

ind = find(zvel(:,2)>-peak_vel_threshold);
zvel(ind,:) = [];
zvel(length(zvel),:)=[];
vel= zvel(:,2);

np_for_envelope=38; % A MODULER
[upper,lo] = envelope(vel,np_for_envelope,'peak');

figure(2);
hold on ;
plot(zvel(:,1),-lo,'r-');
plot(zvel(:,1),-zvel(:,2),'o-');

avg = mean(-lo)
sigma = std(-lo)
percent_uneven=sigma/avg*100
skew= skewness(-lo)
max1=-min(lo)
min1=-max(lo)
%
% a=[1.487 1.952 2.344];% just the peaks
% avg = mean(a)
% sigma = std(a)
% percent_uneven=sigma/avg*100
% skew= skewness(a)
% max1=-min(a)
% min1=-max(a)


grid on
legend('peak envelope','Z vel signal')
xlabel('Distance (m)')
ylabel('Z-Velocity')