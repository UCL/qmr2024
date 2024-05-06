function qmrdiff2024
% QMRDIFF Quantitative MRI Diffusion Challenge 2024
%
% David Atkinson  D.Atkinson@ucl.ac.uk

% Edit for your location of data
dinput = '/Users/davidatkinson/Library/CloudStorage/OneDrive-UniversityCollegeLondon/data/DiffusionChallenge/00005-MR-ep2d_basic_clinical_ORIG' ;

dlist = dir(fullfile(dinput,'*')) ;

ifn = 0 ;
ffns = {} ;

for id = 1:length(dlist)
    if ~dlist(id).isdir && isdicom(fullfile(dlist(id).folder, dlist(id).name))
        ifn=ifn+1;
        ffns{ifn} = fullfile(dlist(id).folder, dlist(id).name) ;
    end
end

dinfo = dfparse(ffns) ;
bvs = unique([dinfo.DiffusionBValue]) ;
bvs(bvs==0) = [] ;

% d2mat returns volume and information about data
[vb0, mb0]    = d2mat(dinfo,{'slice','bv'},'bv',0,'op','fp') ;
[vbnz, mbnz] = d2mat(dinfo,{'slice','bv','bdirec','ddty'},'bv',bvs,'bdirec',[1 2 3],'ddty',1,'op','fp') ;

% Combine the directions using geomteric mean of signals
% vbnz_iso
vbiso = geomean(vbnz,5) ; % [ny nx nslice nbv]

% prepend the b=0 data
vbiso = cat(4,vb0,vbiso) ;
bvs = [0 bvs] ;

% calcADC is in another repo but copied below for info.
[ADC, S0] = calcADC(vbiso, bvs) ;
sviewer(ADC*1e6,mbnz)
sviewer(vbiso,mbnz)

% Load the rois, saved manually as cell array [nslice nroi]
load roipos.mat  % loads pos

hswarm = figure ;
hax = gca ;
hax.YScale = 'log';
hold on
grid on
xlabel('b-value')
ylabel('Log Signal')
title('Log signal against b-value for ROIs in slices 3 and 4')

col{1} = [1 0 0] ;
col{2} = [0.9 0.6 0.13] ;
col{3} = [0.49 0.18 0.56] ;
col{4}= [0.46 0.67 0.19] ;
col{5}= [0 0.45 0.75] ;

hptot = [] ;
bvalx = [0:10:2000] ;
%for islice = 1:size(pos,1)
for islice = [3 4]
    if ~isempty(pos{islice,1})

        UD = sviewer(ADC*1e6,mbnz,"indexD3",islice) ;
        hind = sub2ind([UD.nslice UD.nd4],UD.dslice, UD.dd4) ;
        if isfield(UD,'probe')
            cdata = UD.probe(:, :, UD.dslice, UD.dd4) ;
        else
            cdata = UD.him(hind).CData ;
        end
        for iroi = 1:5
            r_this = drawpolygon(UD.ha, "Position",pos{islice,iroi}) ;
            bw = createMask(r_this, UD.him(hind)) ;
            dat = cdata(bw) ;
            disp("Slice: "+islice+" roi: "+iroi+" Mean: "+mean(dat(:))+...
                " std: "+std(dat(:)))
            madc = mean(dat(:)) ;
            S0slice = S0(:,:,islice) ;
            S0dat = S0slice(bw) ;
            mS0 = mean(S0dat(:)) ;
            spred = mS0*exp(-madc*bvalx*1e-6) ;
            if rem(islice,2) == 0
                lspec = '--' ;
            else
                lspec = ':' ;
            end

            hp = plot(hax,bvalx,spred','Color',col{iroi},'LineWidth',2, ...
                'LineStyle', lspec, ...
                'DisplayName',['Slice: ',num2str(islice),' ROI: ',num2str(iroi), ...
                ' ADC: ',num2str(mean(dat(:)),'%4.1f'),' (',num2str(std(dat(:)),2),')']);

            hptot = [hptot hp] ;

            for ibv = 1:length(bvs)
                bdata = vbiso(:,:,islice,ibv) ;
                yval = bdata(bw) ;
                s= swarmchart(hax, repmat(bvs(ibv),[length(yval(:)) 1]),yval(:),...
                    'XJitter','density','XJitterWidth',100,'CData',col{iroi}, ...
                    'DisplayName','') ;
            end
        end
        close(UD.hf)
    end
end
legend(hptot)


hmap=figure;
hamap = gca ;
hold on
sl=3 ;
imshow(ADC(:,:,sl)*1e6,'XData',mbnz.geom(sl).XData,'YData',mbnz.geom(sl).YData, ...
    'Parent',hamap,'DisplayRange',[0 2000])
for iroi=1:5
  r_this = drawpolygon(hamap, "Position",pos{sl,iroi},'Color',col{iroi}, ...
      'FaceAlpha',0) ;
end

cb = colorbar(hamap, 'southoutside', 'AxisLocation','in', 'FontSize',12) ;
cmap = colorcet('L16') ;
colormap(cmap) ;

% Copy of calcADC used here.
% function [ADC, S0] = calcADC(volb, bvVec, bvUse)
% % calcADC Calculate ADC by linear fitting to log of data (multiple b
% % values)
% %
% % [ADC, S0] = calcADC(volb, bvVec)
% % [ADC, S0] = calcADC(volb, bvVec, bvUse)
% %
% % volb [ny nx nz nbv]  or  [ny nx nbv]
% % bvVec should have length nbv
% % bvUse defaults to bvVec, otherwise only the bvalues present in bvUse will
% % be included in the ADC calculation
% %
% % Solves the following for S0 and ADC using a least squares fit 
% % to the log of the data:
% %   S=S0 exp(-ADC.b)
% % i.e. 
% % ln(S) = ln(S0) - ADC.b
% % ln(S) = ADC.(-b) + ln(S0)
% % Y = mX + c
% % Y = [X 1] [m c]'
% % X = A\B ;
% % ADC = X(1,:) ;
% % S0 = X(2,:) ;
% % 
% % Examples
% % dinfo = datparse ;
% % [volb, matp] = d2mat(dinfo,{'slice','bv'},'op','fp') ;
% % [ADC, S0] = calcADC(volb, matp.bvVec) ;
% %  displayADC(ADC*1e6, 2000)
% %
% % dinfo = dfparse(ffn) ;
% % [vb0, mb0] = d2mat(dinfo,{'slice','bv'},'bv',0,'op','fp%) ;
% % bvs = unique([dinfo.DiffusionBValue]) ;
% % [vbiso, mbiso] = d2mat(dinfo,{'slice','bv','ddty'},'ddty',2,'bv',bvs(2:end),'op','fp%) ;
% % volb = cat(4, vb0, vbiso);
% % [ADC, S0] = calcADC(volb, bvs) ;
% %
% %
% %
% % David Atkinson
% %
% % See also dwi2ADC
% 
% % Reshape so that single slice is still 4D (inc bvaues)
% ndv = ndims(volb) ;
% if ndv==3
%     volb = reshape(volb,[size(volb,1) size(volb,2) 1 ndv]) ;
% end
% 
% volb(~isfinite(volb)) = 0 ;
% loc = volb < 0 ;
% if ~isempty(loc)
%     warning('Negative signal in input to calcADC')
%     volb(loc) = 0 ;
% end
% 
% 
% % Restrict bvalues if requested.
% if nargin == 3
%    [bvIntersect, ia] = intersect(bvVec,bvUse) ;
%    if length(bvUse)~=length(bvIntersect) 
%        warning(['Not all b-values in bvUse are in data'])
%    end
% 
%    bvVec = bvVec(ia) ;
%    volb = volb(:,:,:,ia) ;
% end
% 
% [ny nx nz nbv] = size(volb) ;
% if nbv ~= length(bvVec) 
%     error(['Last dimension of volb must be same as number of bvalues'])
% end
% 
% A = cat(2, -bvVec(:), ones([length(bvVec) 1])) ;
% 
% B = log(volb) ;
% B = reshape(B,[ny*nx*nz nbv]) ;
% 
% Btest = sum(B,2) ;
% [row] = find(~isfinite(Btest)) ;
% B(row, :) = 0 ;
% 
% B = B.' ;
% 
% X = A\B ;
% 
% ADC = X(1,:) ;
% S0 = X(2,:) ;
% 
% ADC(row) = 0 ;
% S0(row) = 0 ;
% 
% ADC = reshape(ADC,[ny nx nz]) ;
% S0 = exp(reshape(S0,[ny nx nz])) ;


