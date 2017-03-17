function [Nd,Coord,Disp,Mass]=node_plot_3Df(fname,fpath)
% Script to plot a 3D model
% Written: Gerard J. O'Reilly
% Last Updated: Feb 2015

% close all; clear all; clc;
% [fname,fpath]=uigetfile('*.txt','Select node coordinates file  ');

n=linecount(fullfile(fpath,fname));
fid=fopen(fullfile(fpath,fname),'r');
out=textscan(fid,'%s','delimiter','\n');
fclose(fid);

Nd=[];
Coord=[];
Disp=[];
Mass=[];
for i=1:n-1
    a=strmatch('Node',out{1,1}(i,:)); % Look for Nodes
    if isempty(a)==0
        % Found an node so extract more info
%         fprintf('found Node at line: %d\n',i);
        temp1=textscan(out{1}{i},' Node: %d');
        Nd=[Nd; temp1{1}];
    end
    b=strmatch('Coordinates',out{1,1}(i,:)); % Look for Coordinates
    if isempty(b)==0  
        temp2=textscan(out{1}{i},'	Coordinates  : %f %f %f ');
        Coord=[Coord; [temp2{1}, temp2{2}, temp2{3}]];
    end
    c=strmatch('Disps',out{1,1}(i,:)); % Look for Disps
    if isempty(c)==0
        temp3=textscan(out{1}{i},'	Disps: %f %f %f %f %f %f \n');
        Disp=[Disp; [temp3{1}, temp3{2}, temp3{3}]];
    end
    d=strmatch('Mass',out{1,1}(i,:)); % Look for Mass
    if isempty(d)==0
        temp4=textscan(out{1}{i+1},'%f 0 0 0 0 0 ');
        temp5=textscan(out{1}{i+2},'0 %f 0 0 0 0 ');
        temp6=textscan(out{1}{i+3},'0 0 %f 0 0 0 ');
        temp7=textscan(out{1}{i+4},'0 0 0 %f 0 0 ');
        temp8=textscan(out{1}{i+5},'0 0 0 0 %f 0 ');
        temp9=textscan(out{1}{i+6},'0 0 0 0 0 %f ');
        Mass=[Mass; [temp4{:}, temp5{:}, temp6{:}, temp7{:}, temp8{:}, temp9{:}]];
    end

end


