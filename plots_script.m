figure

g = 4;

for j = 1:g
    subplot(g/2,2,j)
    contourf(M(:,:,j*floor(length(t)/g)),'ShowText','on','LevelStepMode',...
        'manual','LevelStep',tenCont);
    xlabel('x position (node)'); ylabel('y position (node)')
    lbl = num2str(t(j*floor(length(t)/g)));
    lbl = [lbl,' seconds'];
    ttl = ['Time = ',lbl]; title(ttl);
    set(gca,'Ydir','reverse');
    axis equal
end
pathStr = '/Users/jonathon/Documents/Thesis/GitRepo/Thermal-Model/latex/figures/';
tempInitStr = [',it' num2str(tempInit)];
powerInputStr = [',pi' num2str(Qs)];
nameStr = 'Transient';
pathStr = strcat(pathStr,simTitle,nameStr,tempInitStr,powerInputStr);
if strcmp(prnt,'Y') == 1
    print(pathStr,'-depsc')
else
end

figure

contourf(M(:,:,length(t)-1),10,'ShowText','on');
grid on
xlabel('x position (node)'); ylabel('y position (node)')
lbl = num2str(t(end));
lbl = [lbl,' seconds'];
ttl = ['Temperature Distribution at ',lbl]; title(ttl);
set(gca,'Ydir','reverse');
axis equal
pathStr = '/Users/jonathon/Documents/Thesis/GitRepo/Thermal-Model/latex/figures/';
nameStr = 'Dist';
pathStr = strcat(pathStr,simTitle,nameStr,tempInitStr,powerInputStr);
if strcmp(prnt,'Y') == 1
    print(pathStr,'-depsc')
else
end

figure
plot(t,y(:,1),'LineWidth',2)
hold on
plot(t,y(:,n),'-.','LineWidth',2)
plot(t,y(:,m*n),'LineWidth',2)
plot(t,y(:,m*n-(n-1)),'--','LineWidth',2)
plot(t,y(:,ceil(m*n/2)),'LineWidth',2)
xlabel('Time (seconds)'); ylabel('Temperature (Degrees Celcius)');
legend('Bottom Left','Top Left','Top Right','Bottom Right','Center Element')
title('Temperature Response')
pathStr = '/Users/jonathon/Documents/Thesis/GitRepo/Thermal-Model/latex/figures/';
nameStr = 'Response';
pathStr = strcat(pathStr,simTitle,nameStr,tempInitStr,powerInputStr);
if strcmp(prnt,'Y') == 1
    print(pathStr,'-depsc')
else
end

timeSim/60;
timeSimSec = num2str(timeSim,3);
timeSimMin = num2str(timeSim/60,3);
if timeSim <= 60
    disp(['Simulation Run time = ',timeSimSec,' seconds'])
else
    disp(['Simulation Run time = ',timeSimMin,' minutes'])
end

figure
for i = 0:n-1
    powerMap(:,i+1) = B(i*n+1:i*n+n,1);
end
powerMap = [powerMap zeros(n,1)]; powerMap = [powerMap;zeros(1,n+1)];
pcolor(powerMap);
colormap(gray(2));
colormap(flipud(colormap));
xlabel('x position (node)'); ylabel('y position (node)')
set(gca,'Ydir','reverse'); set(gca,'YTick',[]); set(gca,'XTick',[]);
axis square
title('Power Input Map')
pathStr = '/Users/jonathon/Documents/Thesis/GitRepo/Thermal-Model/latex/figures/';
nameStr = 'Map';
pathStr = strcat(pathStr,simTitle,nameStr,tempInitStr,powerInputStr);
if strcmp(prnt,'Y') == 1
    print(pathStr,'-depsc')
else
end