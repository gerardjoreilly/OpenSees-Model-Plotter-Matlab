function []=model_plot(dims,vw,Npf,Epf,Xamp,Yamp,Zamp,eig)
%% Function to plot either a 2D or 3D model from OpenSees
% Copyright by Gerard J. O'Reilly, 2017
% Written: Gerard J. O'Reilly
% Date: February 2014
% Last Updated: March 2017
%
% ---------------- INPUTS ------------------
% model_plot(dims,vw,Npf,Epf,Xamp,Yamp,Zamp,eig)
% dims: 	no. of dimensions in model (2 or 3)
% vw:       View Plane (xy, xz, yz or 3d)
% Npf:      node label plot flag (1 for labels, 0 for none)
% Epf:      element label plot flag (1 for labels, 0 for none)
% Xamp:     Amplification on Xcoord disps (3D only, omit for 2D) 
% Yamp:     Amplification on Ycoord disps (3D only, omit for 2D)
% Zamp:     Amplification on Zcoord disps (3D only, omit for 2D)
% eig:      number of mode shapes, 0 to do nothing (3D only)

% Requires the complete nodes and elemenst printout from OpenSees
% which is obtained by just writing
% print model.txt


%% Plotting variables
lw1=2.0;
lw2=1.5;
ms1=3;
ms2=5;
fs1=14;

%% Arguments Stuff
if nargin==0
    fprintf('model_plot(dims,vw,Npf,Epf,Xamp,Yamp,Zamp,eig)\n');
    fprintf('dims: 	no. of dimensions in model (2 or 3)\n');
    fprintf('vw:       View (xy, xz, yz, 3d)\n');
    fprintf('Npf:      node label plot flag (1 for labels, 0 for none)\n');
    fprintf('Epf:      element label plot flag (1 for labels, 0 for none)\n');
    fprintf('Xamp:     Amplification on Xcoord disps (3D only, omit for 2D)\n');
    fprintf('Yamp:     Amplification on Ycoord disps (3D only, omit for 2D)\n');
    fprintf('Zamp:     Amplification on Zcoord disps (3D only, omit for 2D)\n');
    fprintf('eig:      number of mode shapes, 0 to do nothing (3D only)\n');
    error('No input args');
end
if nargin<5
    Xamp=1;
    Yamp=1;
    Zamp=1;
    eig=0;
end

%% Prompt for the input files
fprintf('----- 1. Select model file\n');
fprintf('----- 2. Select eigenvector file\n');
fprintf('----- 3. Select periods file\n');

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
        bb=strfind(out{1,1}(i,:),'ForceBeamColumn2d'); % Look for 2d beam columns
        if isempty(bb{1,1})==0
            % found a 2d beam column so extract more info
%             fprintf('found ForceBeamColumn2d at line: %d\n',i);
            temp1=textscan(out{1}{i},'Element: %d Type: ForceBeamColumn2d Connected Nodes: %d %d');
            ele=[ele; temp1{1}];
            iNd=[iNd; temp1{2}];
            jNd=[jNd; temp1{3}];
            Typ{end+1,1}='ForceBeamColumn2d';
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
        cc=strfind(out{1,1}(i,:),'Truss'); % Look for Truss elements
        if isempty(cc{1,1})==0
            % found a Truss so extract more info
%             fprintf('found Truss at line: %d\n',i);
            temp1=textscan(out{1}{i},'Element: %d type: Truss  iNode: %d jNode: %d Area: %f Mass/Length: %f cMass: %f ');
            ele=[ele; temp1{1}];
            iNd=[iNd; temp1{2}];
            jNd=[jNd; temp1{3}];
            Typ{end+1,1}='Truss';
        end
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
    e=strmatch('ElasticBeam2d',out{1,1}(i,:)); % Look for Elastic Elements
    if isempty(e)==0
%         fprintf('found ElasticBeam2d at line: %d\n',i);
        temp5=textscan(out{1}{i},'ElasticBeam2d: %d'); 
        temp6=textscan(out{1}{i+1},'Connected Nodes: %d %d'); 
        ele=[ele; temp5{1}];
        iNd=[iNd; temp6{1}];
        jNd=[jNd; temp6{2}];
        Typ{end+1,1}='ElasticBeam2d';
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
    
end

%% Condense all together
elms=[ele, iNd, jNd];
nelms=length(elms);

%% Get the node coordinates
if dims ==2
    [node, coord, disp, mass]=node_plot_2Df(fname,fpath);
elseif dims ==3
    [node, coord, disp, mass]=node_plot_3Df(fname,fpath); % Updated function Feb 2015
end

if sum(disp)==0
    disp=zeros(length(node),3);
end

%% Get the Eigenmodes
if eig>0
    periods=load(fullfile(fpath_p,fname_p));
    if dims==2
        [~,~,eDisp]=eigennode_plot_2Df(fname_en,fpath_en,eig);
    elseif dims==3
        [~,~,eDisp]=eigennode_plot_3Df(fname_en,fpath_en,eig);
    end
end

%% Plot the model
if dims==2
    h1=figure;
    hold on;
    for i=1:nelms
        if strcmp(Typ{i,1},'ForceBeamColumn2d')==1
            clr='r'; lw=lw1;
        elseif strcmp(Typ{i,1},'ElasticBeam2d')==1
            clr='b'; lw=lw2;
        elseif strcmp(Typ{i,1},'ZeroLength')==1
            clr='y'; lw=lw1;
        elseif strcmp(Typ{i,1},'CorotTrussSection')==1
            clr='g'; lw=lw1;
        elseif strcmp(Typ{i,1},'Truss')==1
            clr='r'; lw=lw1;
        else 
            clr='k'; lw=lw1;
        end
        k1=find(node(:,1)==elms(i,2));
        k2=find(node(:,1)==elms(i,3));
        x1=(coord(k1,1)+disp(k1,1)*Xamp); y1=(coord(k1,2)+disp(k1,2)*Yamp);
        x2=(coord(k2,1)+disp(k2,1)*Xamp); y2=(coord(k2,2)+disp(k2,2)*Yamp);
        xmp=mean([x1, x2]); ymp=mean([y1, y2]);
        plot([x1; x2],[y1; y2],'-','color',clr,'linewidth',lw);
        if Epf == 1
            text(xmp,ymp,num2str(elms(i,1)),'fontsize',fs1,'EdgeColor','k','BackgroundColor','c');
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
        plot(coord(i,1)+disp(i,1)*Xamp,coord(i,2)+disp(i,2)*Yamp,'s','MarkerFaceColor',clrn,'MarkerEdgeColor',clrn,'markersize',msn);
        if Npf == 1
            text(coord(i,1)+disp(i,1)*Xamp+0.1,coord(i,2)+disp(i,2)*Yamp,num2str(node(i,1)),'fontsize',fs1,'EdgeColor','k','BackgroundColor','y');
        end
    end
    xlabel('X-Coordinates');
    ylabel('Y-Coordinates'); 
    if sflag==1
        sPath=uigetdir;
        sName=input('Enter filename: ','s');
        save2pdf(fullfile(sPath,sName),h1);
    end
elseif dims==3
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
    xlabel('X-Coordinates'); xlim([min(coord(:,1)+disp(:,1)*Xamp) max(coord(:,1)+disp(:,1)*Xamp)]);
    ylabel('Y-Coordinates'); 
    zlabel('Z-Coordinates'); zlim([min(coord(:,3)+disp(:,3)*Zamp) max(coord(:,3)+disp(:,3)*Zamp)]);
end

%% Eigenmodes
if eig>0
    if dims ==2
        for mds=1:eig
            figure; hold on; box on;
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
                x1=(coord(k1,1)+eDisp{mds,1}(k1,1)*Xamp); y1=(coord(k1,2)+eDisp{mds,1}(k1,2)*Yamp); 
                x2=(coord(k2,1)+eDisp{mds,1}(k2,1)*Xamp); y2=(coord(k2,2)+eDisp{mds,1}(k2,2)*Yamp);
                xmp=mean([x1, x2]); ymp=mean([y1, y2]);
                plot([x1; x2],[y1; y2],'-','color',clr,'linewidth',lw);
                if Epf == 1
                    text(xmp,ymp,num2str(elms(i,1)),'fontsize',fs1,'EdgeColor','k','BackgroundColor','c');
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
                plot(coord(i,1)+eDisp{mds,1}(i,1)*Xamp,coord(i,2)+eDisp{mds,1}(i,2)*Yamp,'s','MarkerFaceColor',clrn,'MarkerEdgeColor',clrn,'markersize',msn);
                if Npf == 1
                    text(coord(i,1)+eDisp{mds,1}(i,1)*Xamp+0.1,coord(i,2)+eDisp{mds,1}(i,2)*Yamp,num2str(node(i,1)),'fontsize',fs1,'EdgeColor','k','BackgroundColor','y');
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
            xlabel('X'); %xlim([min(coord(:,1)+eDisp{mds,1}(:,1)*Xamp) max(coord(:,1)+eDisp{mds,1}(:,1)*Xamp)]);
            ylabel('Y'); 
            title(sprintf('Mode %d  T=%1.3fs',mds,periods(mds)),'fontsize',20);
        end
    elseif dims==3
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
                x1=(coord(k1,1)+eDisp{mds,1}(k1,1)*Xamp); y1=(coord(k1,2)+eDisp{mds,1}(k1,2)*Yamp); z1=(coord(k1,3)+eDisp{mds,1}(k1,3)*Zamp); 
                x2=(coord(k2,1)+eDisp{mds,1}(k2,1)*Xamp); y2=(coord(k2,2)+eDisp{mds,1}(k2,2)*Yamp); z2=(coord(k2,3)+eDisp{mds,1}(k2,3)*Zamp);
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
                plot3(coord(i,1)+eDisp{mds,1}(i,1)*Xamp,coord(i,2)+eDisp{mds,1}(i,2)*Yamp,coord(i,3)+eDisp{mds,1}(i,3)*Zamp,'s','MarkerFaceColor',clrn,'MarkerEdgeColor',clrn,'markersize',msn);
                if Npf == 1
                    text(coord(i,1)+eDisp{mds,1}(i,1)*Xamp+0.1,coord(i,2)+eDisp{mds,1}(i,2)*Yamp,coord(i,3)+eDisp{mds,1}(i,3)*Zamp+0.1,num2str(node(i,1)),'fontsize',fs1,'EdgeColor','k','BackgroundColor','y');
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
            xlabel('X-Coordinates'); %xlim([min(coord(:,1)+eDisp{mds,1}(:,1)*Xamp) max(coord(:,1)+eDisp{mds,1}(:,1)*Xamp)]);
            ylabel('Y-Coordinates'); 
            zlabel('Z-Coordinates'); %zlim([min(coord(:,3)+eDisp{mds,1}(:,3)*Zamp) max(coord(:,3)+eDisp{mds,1}(:,3)*Zamp)]);
            title(sprintf('Mode %d  T=%1.3fs f=%1.3fs',mds,periods(mds),1/periods(mds)),'fontsize',20);
        end
    end
end











