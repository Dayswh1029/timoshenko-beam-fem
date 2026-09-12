
k_arr=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40];
l=length(k_arr);

j_arr=['m', 'd', 'sn'];
j_col=['g','k','r'];

for k=1:l
filename= strcat('Sq',num2str(k_arr(k)),'_');
figure(k)
for j=1:3
if j==3
    filename1= strcat(filename,'sn','.mat');
else
filename1= strcat(filename,j_arr(j),'.mat');
end
load(filename1);


x=(FRF_Point2_Point1.x_values.start_value:FRF_Point2_Point1.x_values.increment:(FRF_Point2_Point1.x_values.start_value+FRF_Point2_Point1.x_values.increment*(FRF_Point2_Point1.x_values.number_of_values-1)));
y=abs(FRF_Point2_Point1.y_values.values);

rr=strcat(j_col(j),'--');
loglog(x,y,rr);
hold on; grid on;

end
legend ('FRF m','FRF d', 'FRF sn');
title(filename);
end

