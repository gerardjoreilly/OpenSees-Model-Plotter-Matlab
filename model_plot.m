function [T,phi_mds,Mnstar]=model_plot(vw,Npf,Epf,Xamp,Yamp,Zamp,eig,node_list)
%% Function to plot a 3D model from OpenSees
% Copyright by Gerard J. O'Reilly, 2017
% Written: Gerard J. O'Reilly
% Date Created: February 2014
%
% ---------------- INPUTS ------------------
% [T,phi_mds,Mnstar]=model_plot(vw,Npf,Epf,Xamp,Yamp,Zamp,node_list)
% vw:       View Plane (xy, xz, yz or 3d)
% Npf:      node label plot flag (1 for labels, 0 for none)
% Epf:      element label plot flag (1 for labels, 0 for none)
% Xamp:     Amplification on Xcoord disps
% Yamp:     Amplification on Ycoord disps
% Zamp:     Amplification on Zcoord disps
% eig:      number of mode shapes, 0 to do nothing
% node_list:list of nodes to return the modal displacements at, returns a
%           cell with eig number of cells and node_list x ndf array in each

% clear;
% clc;
% vw='3d';
% Npf=0;
% Epf=0;
% Xamp=1;
% Yamp=1;
% Zamp=1;
% eig=0;
% mds=[];
% nargin=10;

%% Plotting variables
lw1=2.0;
lw2=1.5;
ms1=3;
ms2=5;
fs1=14;

%% Arguments Stuff
if nargin==0
    fprintf('[T,phi_mds,Mnstar]=model_plot(vw,Npf,Epf,Xamp,Yamp,Zamp,eig,node_list)\n');
    fprintf('vw:       View (xy, xz, yz, 3d)\n');
    fprintf('Npf:      node label plot flag (1 for labels, 0 for none)\n');
    fprintf('Epf:      element label plot flag (1 for labels, 0 for none)\n');
    fprintf('Xamp:     Amplification on Xcoord disps\n');
    fprintf('Yamp:     Amplification on Ycoord disps\n');
    fprintf('Zamp:     Amplification on Zcoord disps\n');
    fprintf('eig:      number of mode shapes, 0 to do nothing\n');
    fprintf('node_list:list of nodes to return the modal displacements at, returns a\n');
    fprintf('          cell with eig number of cells and node_list x ndf array in each\n\n');
    
    error('No input args');
end
if nargin<4
    Xamp=1;
    Yamp=1;
    Zamp=1;
    eig=0;
    mds=[];
end

%% Prompt for the input files
fprintf('----- 1. Select model file\n');
if eig>0
    fprintf('----- 2. Select eigenvector file\n');
    fprintf('----- 3. Select periods file\n');
end

[fname,fpath]=uigetfile('*.txt','Select model:');

if eig>0
    [fname_en,fpath_en]=uigetfile('*.txt','Select eigen file  ');
    [fname_p,fpath_p]=uigetfile('*.txt','Select periods file  ');
end

%% Open the elements file and get stuff
n=linecount(fullfile(fpath,fname));
fid=fopen(fullfile(fpath,fname),'r');
out=textscan(fid,'%s','delimiter','\n');
fclose(fid);

ele=[];
iNd=[];
jNd=[];
Typ={};
for i=1:n-1
    a=strmatch('Element',out{1,1}(i,:)); % Look for Elements
    if isempty(a)==0
        % Found an element so extract more info
        b=strfind(out{1,1}(i,:),'ForceBeamColumn3d'); % Look for 3d beam columns
        if isempty(b{1,1})==0
            % found a 3d beam column so extract more info
%             fprintf('found ForceBeamColumn3d at line: %d\n',i);
            temp1=textscan(out{1}{i},'Element: %d Type: ForceBeamColumn3d Connected Nodes: %d %d');
            ele=[ele; temp1{1}];
            iNd=[iNd; temp1{2}];
            jNd=[jNd; temp1{3}];
            Typ{end+1,1}='ForceBeamColumn3d';
        end

        c=strfind(out{1,1}(i,:),'ZeroLength'); % look for zerolengths 
        if isempty(c{1,1})==0
            % found a zerolength so extract
%             fprintf('found ZeroLength at line: %d\n',i);
            temp2=textscan(out{1}{i},'Element: %d type: ZeroLength iNode: %d jNode: %d');
            ele=[ele; temp2{1}];
            iNd=[iNd; temp2{2}];
            jNd=[jNd; temp2{3}]; 
            Typ{end+1,1}='ZeroLength';
        end
%         cc=strfind(out{1,1}(i,:),'Truss'); % Look for Truss elements
%         if isempty(cc{1,1})==0
%             % found a Truss so extract more info
% %             fprintf('found Truss at line: %d\n',i);
%             temp1=textscan(out{1}{i},'Element: %d type: Truss  iNode: %d jNode: %d Area: %f Mass/Length: %f cMass: %f ');
%             ele=[ele; temp1{1}];
%             iNd=[iNd; temp1{2}];
%             jNd=[jNd; temp1{3}];
%             Typ{end+1,1}='Truss';
%         end
    end
    d=strmatch('ElasticBeam3d',out{1,1}(i,:)); % Look for Elastic Elements
    if isempty(d)==0
%         fprintf('found ElasticBeam3d at line: %d\n',i);
        temp3=textscan(out{1}{i},'ElasticBeam3d: %d'); 
        temp4=textscan(out{1}{i+1},'Connected Nodes: %d %d'); 
        ele=[ele; temp3{1}];
        iNd=[iNd; temp4{1}];
        jNd=[jNd; temp4{2}];
        Typ{end+1,1}='ElasticBeam3d';
    end

    f=strmatch('CorotTrussSection',out{1,1}(i,:)); % Look for CorotTrussSection Elements
    if isempty(f)==0
%         fprintf('CorotTrussSection at line: %d\n',i);
        temp5=textscan(out{1}{i},'CorotTrussSection, tag: %d'); 
        temp6=textscan(out{1}{i+1},'Connected Nodes: %d %d'); 
        ele=[ele; temp5{1}];
        iNd=[iNd; temp6{1}];
        jNd=[jNd; temp6{2}];
        Typ{end+1,1}='CorotTrussSection';
    end
    
    g=strfind(out{1,1}(i,:),'TrussSection'); % Look for TrussSection Elements
    if g{:}>0
        fprintf('TrussSection at line: %d\n',i);
        temp5=textscan(out{1}{i},'Element: %d type: TrussSection  iNode: %d jNode: %d Mass density/length: %f cMass: %f ');  
        ele=[ele; temp5{1}];
        iNd=[iNd; temp5{2}];
        jNd=[jNd; temp5{3}];
        Typ{end+1,1}='TrussSection';
    end
    
end

%% Condense all together
elms=[ele, iNd, jNd];
nelms=length(elms);

%% Get the node coordinates
[node, coord, disp, mass]=node_plot_3Df(fname,fpath);

% Needed this for debugging
if sum(disp)==0
%     warning('Model building/debugging mode');
    disp=zeros(length(node),3);
end

%% Get the modal properties
% This will need the output of the modalAnalysis.tcl procedure
if eig>0
    % Get the periods
    T=load(fullfile(fpath_p,fname_p));
    
    % Get the mode shapes
    [~,~,phi]=eigennode_plot_3Df(fname_en,fpath_en,eig);
end

%% Extract the mode shapes
if eig>0 && isempty(node_list)==0 % Do if requested
   phi_mds=cell(eig,1); % Initialise the cell
   for i=1:eig
       for j=1:length(node_list)
           phi_mds{i}(j,:)=phi{i}(find(node_list(j)==node),:);
       end
   end
end

%% Do some computation
% Set up the mass matrix
[lmass,~]=size(mass);
M=diag(reshape(mass',lmass*6,1));

% Compute the total mass
Mtotal=sum(reshape(mass',lmass*6,1));

% Get the influence vector
temp=mass./mass; % get the dofs with mass
temp(isnan(temp))=0; % remove nans
r=reshape(temp',lmass*6,1);

for mds=1:eig
    % Compute the modal mass
    Mn(mds)=reshape(phi{mds}',lmass*6,1)'*M*reshape(phi{mds}',lmass*6,1);

    % Compute modal participation factor
    Ln(mds)=reshape(phi{mds}',lmass*6,1)'*M*r;
    
    % Compute the normalised effective modal mass
    Mnstar(mds)=Ln(mds)^2/Mn(mds)/Mtotal;
end

%% Plot the model
h2=figure;
hold on;
for i=1:nelms
    if strcmp(Typ{i,1},'ForceBeamColumn3d')==1
        clr='b'; lw=lw1;
    elseif strcmp(Typ{i,1},'ElasticBeam3d')==1
        clr='k'; lw=lw2;
    elseif strcmp(Typ{i,1},'ZeroLength')==1
        clr='r'; lw=lw1;
    elseif strcmp(Typ{i,1},'CorotTrussSection')==1
        clr='g'; lw=lw1;
    elseif strcmp(Typ{i,1},'Truss')==1
        clr=[1 .5 0]; lw=lw1;
    else 
        clr='k'; lw=lw1;
    end

    k1=find(node(:,1)==elms(i,2));
    k2=find(node(:,1)==elms(i,3));
    x1=(coord(k1,1)+disp(k1,1)*Xamp); y1=(coord(k1,2)+disp(k1,2)*Yamp); z1=(coord(k1,3)+disp(k1,3)*Zamp); 
    x2=(coord(k2,1)+disp(k2,1)*Xamp); y2=(coord(k2,2)+disp(k2,2)*Yamp); z2=(coord(k2,3)+disp(k2,3)*Zamp);
    xmp=mean([x1, x2]); ymp=mean([y1, y2]); zmp=mean([z1, z2]);
    plot3([x1; x2],[y1; y2],[z1; z2],'-','color',clr,'linewidth',lw);
    if Epf == 1
        text(xmp,ymp,zmp,num2str(elms(i,1)),'fontsize',fs1,'EdgeColor','k','BackgroundColor','c');
    end
end
for i=1:length(node)
    if isempty(mass)==1
        clrn='k';
        msn=ms1;
    else
        if sum(mass(i,:))>0
            clrn='r';
            msn=ms2;
        else
            clrn='k';
            msn=ms1;
        end
    end
    plot3(coord(i,1)+disp(i,1)*Xamp,coord(i,2)+disp(i,2)*Yamp,coord(i,3)+disp(i,3)*Zamp,'s','MarkerFaceColor',clrn,'MarkerEdgeColor',clrn,'markersize',msn);
    if Npf == 1
        text(coord(i,1)+disp(i,1)*Xamp+0.1,coord(i,2)+disp(i,2)*Yamp,coord(i,3)+disp(i,3)*Zamp+0.1,num2str(node(i,1)),'fontsize',fs1,'EdgeColor','k','BackgroundColor','y');
    end
end
axis equal
if vw=='xy'
    view(0,90);
elseif vw=='xz'
    view(0,0);
elseif vw=='yz'
    view(90,0);
elseif vw=='3d'
    view(-37.5,30);
end
xlabel('X-Coordinates'); xlim([min(coord(:,1)+disp(:,1)*Xamp)-0.0001 max(coord(:,1)+disp(:,1)*Xamp)+0.0001]);
ylabel('Y-Coordinates'); 
zlabel('Z-Coordinates'); zlim([min(coord(:,3)+disp(:,3)*Zamp)-0.0001 max(coord(:,3)+disp(:,3)*Zamp)+0.0001]);

%% Plot the Mode Shapes
if eig>0
    for mds=1:eig
        figure; hold on; box on; grid on;
        for i=1:nelms
            if strcmp(Typ{i,1},'ForceBeamColumn3d')==1
                clr='b'; lw=lw1;
            elseif strcmp(Typ{i,1},'ElasticBeam3d')==1
                clr='k'; lw=lw2;
            elseif strcmp(Typ{i,1},'ZeroLength')==1
                clr='r'; lw=lw1;
            elseif strcmp(Typ{i,1},'CorotTrussSection')==1
                clr='g'; lw=lw1;
            else 
                clr='k'; lw=lw1;
            end
            k1=find(node(:,1)==elms(i,2));
            k2=find(node(:,1)==elms(i,3));
            x1=(coord(k1,1)+phi{mds,1}(k1,1)*Xamp); y1=(coord(k1,2)+phi{mds,1}(k1,2)*Yamp); z1=(coord(k1,3)+phi{mds,1}(k1,3)*Zamp); 
            x2=(coord(k2,1)+phi{mds,1}(k2,1)*Xamp); y2=(coord(k2,2)+phi{mds,1}(k2,2)*Yamp); z2=(coord(k2,3)+phi{mds,1}(k2,3)*Zamp);
            xmp=mean([x1, x2]); ymp=mean([y1, y2]); zmp=mean([z1, z2]);
            plot3([x1; x2],[y1; y2],[z1; z2],'-','color',clr,'linewidth',lw);
            if Epf == 1
                text(xmp,ymp,zmp,num2str(elms(i,1)),'fontsize',fs1,'EdgeColor','k','BackgroundColor','c');
            end
        end
        for i=1:length(node)
            if isempty(mass)==1
                clrn='k';
                msn=ms1;
            else
                if sum(mass(i,:))>0
                    clrn='r';
                    msn=ms2;
                else
                    clrn='k';
                    msn=ms1;
                end
            end
            plot3(coord(i,1)+phi{mds,1}(i,1)*Xamp,coord(i,2)+phi{mds,1}(i,2)*Yamp,coord(i,3)+phi{mds,1}(i,3)*Zamp,'s','MarkerFaceColor',clrn,'MarkerEdgeColor',clrn,'markersize',msn);
            if Npf == 1
                text(coord(i,1)+phi{mds,1}(i,1)*Xamp+0.1,coord(i,2)+phi{mds,1}(i,2)*Yamp,coord(i,3)+phi{mds,1}(i,3)*Zamp+0.1,num2str(node(i,1)),'fontsize',fs1,'EdgeColor','k','BackgroundColor','y');
            end
        end
        axis equal
        if vw=='xy'
            view(0,90);
        elseif vw=='xz'
            view(0,0);
        elseif vw=='yz'
            view(90,0);
        elseif vw=='3d'
            view(-37.5,30);
        end
        xlabel('X-Coordinates'); xlim([min(coord(:,1)+phi{mds,1}(:,1)*Xamp)-0.0001 max(coord(:,1)+phi{mds,1}(:,1)*Xamp)+0.0001]);
        ylabel('Y-Coordinates'); ylim([min(coord(:,2)+phi{mds,1}(:,2)*Yamp)-0.0001 max(coord(:,2)+phi{mds,1}(:,2)*Yamp)+0.0001]);
        zlabel('Z-Coordinates'); zlim([min(coord(:,3)+phi{mds,1}(:,3)*Zamp)-0.0001 max(coord(:,3)+phi{mds,1}(:,3)*Zamp)+0.0001]);
        title(sprintf('Mode %d   T=%1.3fs  f=%1.3fsHz  %%M=%.1f',mds,T(mds),1/T(mds),Mnstar(mds)*100),'fontsize',20);
    end
end

%% Display some info on the command window
if eig>0
    fprintf('%s\n','Mode   T [s]   f [Hz]    %M   Sum %M');
    for mds=1:eig
        fprintf('%s%d%s%.3f%s%.3f%s%.1f%s%.1f\n','  ',mds,'    ',T(mds),'   ',1/T(mds),'    ',Mnstar(mds)*100,'   ',sum(Mnstar(1:mds))*100)
    end
    if isempty(node_list)==0
        for mds=1:eig
            fprintf('%s%d%s\n','Mode: ',mds,' ==============================');
            fprintf('%s\n','DOF:     1      2      3      4      5      6');
            for j=1:length(node_list)
                fprintf('%s%s\n','       ',num2str(phi_mds{mds}(j,:),'%.3f  '));
            end
        end
    end
    
end
