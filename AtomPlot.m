%I used the following code to guess at what is the minimum value. I changed
%the h value every time giving me different results depending on  the h. As
%the guesser it took me about 20 tries to get the potential energy to go down to
%-2.1448 which I figured was a minimum. The h used was 1.6105 and the x,y, and
%z values used were
%([0,0,-.8053,.9953,2.1082],[0,0,0,0,0],[0,2.4782,.5850,0,1.5317])
a0=[0,0,0];
a(:,1)=a0;
vertices=360/(5):360/(5):360;
h=1.6105;
for k=1:5
    a(1,k)=a(1,1)+h*cosd(vertices(k));
    a(3,k)=a(3,1)+h*sind(vertices(k));
end
min=PotentialEnergy(a(1,:),a(2,:),a(3,:)); 
disp(min);
%This is the section for plotting my found minimum using plot3
%plot3([0,0,-.8053,.9953,2.1082],[0,0,0,0,0],[0,2.4782,.5850,0,1.5317]);
plot3([0,0],[0,0],[0,2.4782]);
hold on;
plot3([0,-.8053],[0,0],[0,.5850]);
plot3([0,.9953],[0,0],[0,0]);
plot3([0,2.1082],[0,0],[0,1.5317]);
plot3([0,-.8053],[0,0],[2.4782,.5850]);
plot3([0,.9953],[0,0],[2.4782,0]);
plot3([0,2.1082],[0,0],[2.4782,1.5317]);
plot3([-.8053,.9953],[0,0],[.5850,0]);
plot3([-.8053,2.1082],[0,0],[.5850,1.5317]);
plot3([.9953,2.1082],[0,0],[0,1.5317]);