function [] = pa()
clear
clc
close all
warning off

sim.nx = 100;
sim.ny = 100;
sim.iter = 1000;
sim.dx = 1;
sim.delay = 0.001;
sim.movie = false;
sim.useimboxfilt = true;

board = intializeBoard(sim);
board.set(:,1) = true;
%board.set(:,sim.nx) = true;
board.set(1,:) = true;
board.V(:,1) = 1;

runSimulation(board,sim)
end

function [board] = calcE(board,sim)
board.Ex = zeros(sim.ny-1,sim.nx-1);
board.Ey = zeros(sim.ny-1,sim.nx-1);
for n = 1:(sim.ny-1)
    for m = 1:(sim.nx-1)
        board.Ex(n,m) = -(board.V(n+1,m+1) + board.V(n,m+1) - board.V(n+1,m) - board.V(n,m))/(2*sim.dx);
        board.Ey(n,m) = -(board.V(n+1,m+1) + board.V(n+1,m) - board.V(n,m+1) - board.V(n,m))/(2*sim.dx);
    end
end
board.EMag = sqrt(board.Ex.^2 + board.Ey.^2);
end

function [board] = intializeBoard(sim)
board.V = zeros(sim.ny,sim.nx);
board.set = false(sim.ny,sim.nx);
end

function [] = runSimulation(board,sim)
board = calcE(board,sim);
plt = setupPlot(board,sim);
for n = 1:sim.iter
    if sim.useimboxfilt
        board = performIteration2(board);
    else
        board = performIteration(board);
    end
    if sim.movie
        board = updatePlot(board,plt,sim);
        pause(sim.delay);
    end
end
if ~sim.movie
    updatePlot(board,plt,sim);
end
end

function [board] = updatePlot(board,plt,sim)
board = calcE(board,sim);
set(plt.plotV,'ZData',board.V);
set(plt.plotE,'UData',board.Ex,'VData',board.Ey);
set(plt.plotEMag,'ZData',board.EMag);
set(plt.plotE2,'UData',board.Ex,'VData',board.Ey);
end

function [plt] = setupPlot(board,sim)
x = (0:(sim.nx-1))*sim.dx;
y = (0:(sim.ny-1))*sim.dx;
[X,Y] = meshgrid(x,y);

Ex = ((0:(sim.nx-2))+0.5)*sim.dx;
Ey = ((0:(sim.ny-2))+0.5)*sim.dx;
[EX,EY] = meshgrid(Ex,Ey);

plt.fig = figure('units','normalized','outerposition',[0 0 1 1]);

plt.ax = subplot(1, 2, 1, 'Parent', plt.fig);
hold(plt.ax, 'on');
[~,plt.plotV] = contour(plt.ax,X,Y,board.V,20);
colorbar(plt.ax);
plt.plotE = quiver(plt.ax,EX,EY,board.Ex,board.Ey,10);
xlabel(plt.ax,'x');
ylabel(plt.ax,'y');
title(plt.ax,'Voltage with Electric Field Vectors');
grid(plt.ax,'on');

plt.ax2 = subplot(1, 2, 2, 'Parent', plt.fig);
hold(plt.ax2, 'on');
[~,plt.plotEMag] = contour(plt.ax2,EX,EY,board.EMag,20);
colorbar(plt.ax2);
plt.plotE2 = quiver(plt.ax2,EX,EY,board.Ex,board.Ey,10);
xlabel(plt.ax2,'x');
ylabel(plt.ax2,'y');
title(plt.ax2,'Electric Field Magnitude with Electric Field Vectors');
grid(plt.ax2,'on');
end

function [board] = performIteration2(board)
newV = imboxfilt(board.V,3);
board.V(~board.set) = newV(~board.set);
end

function [board] = performIteration(board)
newV = zeros(size(board.V));
for n = 1:size(board.V,1)
    for m = 1:size(board.V,2)
        if n == 1 || n == size(board.V,1) || m == 1 || m == size(board.V,2)
            if board.set(n,m)
                newV(n,m) = board.V(n,m);
            else
                if n == 1
                    newV(n,m) = board.V(n+1,m);
                elseif n == size(board.V,1)
                    newV(n,m) = board.V(n-1,m);
                elseif m == 1
                    newV(n,m) = board.V(n,m+1);
                elseif m == size(board.V,2)
                    newV(n,m) = board.V(n,m-1);
                end
            end
        else
            if board.set(n,m)
                newV(n,m) = board.V(n,m);
            else
                newV(n,m) = (board.V(n-1,m)+board.V(n+1,m)+board.V(n,m-1)+board.V(n,m+1))/4;
            end
        end
    end
end
board.V = newV;
end