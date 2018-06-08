clc;
close all;
clear all;
home;
opengl hardware

flagPlot = 0; % Boolean equal to 0 for calculations + plots and 1 when only for plots
flagCrop = 1; % To do for Z.ang
loop = struct();

loop(1).pname = 'N:\Projects\2016_Invar_Nutal\2017-10_MTEX\data_min\invar_ascast';
loop(2).pname = 'N:\Projects\2016_Invar_Nutal\2017-10_MTEX\data_min\invar_ascast_tth_full';
loop(3).pname = 'N:\Projects\2016_Invar_Nutal\2017-10_MTEX\data_min\invar_slm_xy';
loop(4).pname = 'N:\Projects\2016_Invar_Nutal\2017-10_MTEX\data_min\invar_slm_xy_tth_relax';
loop(5).pname = 'N:\Projects\2016_Invar_Nutal\2017-10_MTEX\data_min\invar_slm_xy_tth_full';
loop(6).pname = 'N:\Projects\2016_Invar_Nutal\2017-10_MTEX\data_min\invar_slm_z';
loop(7).pname = 'N:\Projects\2016_Invar_Nutal\2017-10_MTEX\data_min\invar_slm_z_tth_relax';
loop(8).pname = 'N:\Projects\2016_Invar_Nutal\2017-10_MTEX\data_min\invar_slm_z_tth_full';
loop(9).pname = 'N:\Projects\2016_Invar_Nutal\2017-10_MTEX\data_min\invar_slm_z_tth_full_2';
loop(10).pname = 'N:\Projects\2016_Invar_Nutal\2017-10_MTEX\data_min\invar_slm_z_tth_homo';
loop(11).pname = 'N:\Projects\2016_Invar_Nutal\2017-10_MTEX\data_min\invar_slm_z_tth_relax_2';

% which files to be imported
loop(1).fname = '2367a.ang';
loop(2).fname = '2365c.ang';
loop(3).fname = 'Y.ang';
loop(4).fname = 'AMY2.ang';
loop(5).fname = 'AMY3.ang';
loop(6).fname = 'Z.ang';
loop(7).fname = 'AMZ2.ang';
loop(8).fname = 'AMZ3.ang';
loop(9).fname = 'VT1-100x.ang';
loop(10).fname = 'VT2-100x.ang';
loop(11).fname = 'VT3-100x.ang';

for jj = 1:11
    clc;
    close all;
    clear color;
    clear CS;
    clear d;
    clear ebsd;
    clear ebsd_smoothed;
    clear ext;
    clear F;
    clear filepath;
    clear fname;
    clear gB;
    clear grain_id;
    clear grains;
    clear grainsCleaned;
    clear grainsStat;
    clear h;
    clear id;
    clear Indexed;
    clear kam_0;
    clear kam_1;
    clear kam_2;
    clear kam_3;
    clear mAngle;
    clear material;
    clear name;
    clear notIndexed;
    clear odf;
    clear ori;
    clear outerBoundary_id;
    clear pfname;
    clear pname;
    clear psi;
    clear psiMean;
    clear tolVal;
    clear tolGB;
    clear toRemove;
    display(jj);
    
    %% Specify Path/File Names
    pname = loop(jj).pname;
    fname = loop(jj).fname;
    pfname = [pname, '\', fname];
    [filepath,name,ext] = fileparts(fname);
    
    if flagPlot
        try
            resMat = load(fullfile([pname,'\results.mat']),'-mat');
            res = resMat.res;
            flagCalc = 1;
        catch
            flagCalc = 0;
        end
    else
        flagCalc = 1;
        res = struct();
        
        %% Define CS and MTEXpref
        % crystal symmetry
        CS = {...
            'notIndexed',...
            crystalSymmetry('432', [2.87 2.87 2.87], 'mineral', 'Iron - Alpha', 'color', 'light blue'),...
            crystalSymmetry('432', [3.56 3.56 3.56], 'mineral', 'Nickel', 'color', 'light green')};
        
        % plotting convention
        setMTEXpref('xAxisDirection','south');
        setMTEXpref('zAxisDirection','outOfPlane');
        
        %% Import the Data
        % create an EBSD variable containing the data
        ebsd = loadEBSD(pfname,CS,'interface','ang',...
            'convertSpatial2EulerReferenceFrame');
        
        %% Grains
        material = 'Nickel';
        ebsd = ebsd(material);
        [grains,ebsd.grainId] = calcGrains(ebsd,'angle',10*degree);
        
        Indexed = grains('Indexed');
        notIndexed = grains('notIndexed');
        
        % Remove not indexed pixels
        if length(notIndexed) > 1
            % the "not indexed grains" we want to remove
            toRemove = notIndexed(notIndexed.grainSize ./ notIndexed.boundarySize<0.8);
            
            % now we remove the corresponding EBSD measurements
            ebsd(toRemove) = [];
            
            % and perform grain reconstruction with the reduces EBSD data set
            %[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd);
            [grains] = calcGrains(ebsd);
        end
        
        % Remove small pixels
        toRemove = Indexed(Indexed.grainSize ./ Indexed.boundarySize<0.8);
        ebsd(toRemove) = [];
        
        if flagCrop
            % Crop ebsd map
            ebsdRaw = ebsd;
            region = [0 0 400 400]; % in micron
            %region = [0 0 800 800]; % in micron
            hold on
            rectangle('position',region,'edgecolor','r','linewidth',2);
            condition = inpolygon(ebsd,region);
            ebsd = ebsd(condition);
        end
        
        tolVal = [1.5 5];
        tolGB = [5 10];
    end
    
    if flagCalc
        for ii = 1:length(tolVal)
            if flagPlot
                CS = res(ii).CS;
                d = res(ii).d;
                ebsd = res(ii).ebsd;
                ebsd_smoothed = res(ii).ebsd_smoothed;
                ext = res(ii).ext;
                F = res(ii).F;
                filepath = res(ii).filepath;
                fname = res(ii).fname;
                gB = res(ii).gB;
                grain_id = res(ii).grain_id;
                grains = res(ii).grains;
                grainsCleaned = res(ii).grainsCleaned;
                grainsStat = res(ii).grainsStat;
                h = res(ii).h;
                id = res(ii).id;
                Indexed = res(ii).Indexed;
                kam_0 = res(ii).kam_0;
                kam_1 = res(ii).kam_1;
                kam_2 = res(ii).kam_2;
                kam_3 = res(ii).kam_3;
                mAngle = res(ii).mAngle;
                material = res(ii).material;
                mdf_mat_mat = res(ii).mdf_mat_mat;
                name = res(ii).name;
                notIndexed = res(ii).notIndexed;
                odf = res(ii).odf;
                odf_mat = res(ii).odf_mat;
                odf_model = res(ii).odf_model;
                odfModel_odf = res(ii).odfModel_odf;
                odf_zero = res(ii).odf_zero;
                ori = res(ii).ori;
                outerBoundary_id = res(ii).outerBoundary_id;
                pf = res(ii).pf;
                pfname = res(ii).pfname;
                pname = res(ii).pname;
                psi = res(ii).psi;
                psiMean = res(ii).psiMean;
                tolVal = res(ii).tol.tolVal;
                tolGB = res(ii).tol.tolGB;
                toRemove = res(ii).toRemove;
            end
            
            if ~flagPlot
                [grains,ebsd.grainId] = calcGrains(ebsd,'angle',tolVal(ii)*degree);
            end
            %         figure; set(gcf, 'Color', [1 1 1]);
            %         plot(grains);
            %         saveFig(pfname, pname, name, [num2str(ii), '_grains_phase_', num2str(tolVal(ii)), 'deg']);
            %
            %% Grain boundaries
            if ~flagPlot
                %gB = grains.boundary;
                gB = grains.boundary(material, material);
                figure; set(gcf, 'Color', [1 1 1]);
            end
            %         plot(grains,'translucent',.3,'micronbar','off');
            %         legend off
            %         hold on
            %         plot(gB,gB.misorientation.angle./degree,'linewidth',1.5);
            %         hold off
            %         mtexColorbar;
            %         %set(mtexColorbar,'Title','misorientation angle');
            %         saveFig(pfname, pname, name, [num2str(ii), '_gb_misor_', num2str(tolVal(ii)), 'deg']);
            %
            if ~flagPlot
                mAngle = gB.misorientation.angle./ degree;
                [~,id] = histc(mAngle,0:tolGB(ii):120);
            end
            
            figure; set(gcf, 'Color', [1 1 1]);
            plot(gB,'linecolor','k');
            hold on
            plot(gB(id==1),'linecolor','b','linewidth',2,'DisplayName',...
                ['<', num2str(tolGB(ii)),'^\circ']);
            hold on;
            plot(gB(id>=2),'linecolor','r','linewidth',2,'DisplayName',...
                ['>', num2str(tolGB(ii)),'^\circ']);
            saveFig(pfname, pname, name, [num2str(ii),'_gb_misor_map_', num2str(tolVal(ii)), 'deg']);
            
            figure; set(gcf, 'Color', [1 1 1]);
            plotAngleDistribution(gB.misorientation);
            xlim([0 120]);
            saveFig(pfname, pname, name, [num2str(ii),'_gb_misor_dist_', num2str(tolVal(ii)), 'deg']);
            
            %         figure; set(gcf, 'Color', [1 1 1]);
            %         plotAngleDistribution(gB.misorientation,'DisplayName','From EBSD');
            %         xlim([0 120]); hold on;
            %         legend('-dynamicLegend','Location','northwest') % update legend
            %
            %         % compute the ODF
            %         if ~flagPlot
            %             odf_mat = calcODF(ebsd(material).orientations,'Fourier');
            %             % compute the uncorrelated material to material MDF
            %             mdf_mat_mat = calcMDF(odf_mat,odf_mat);
            %         end
            %         plotAngleDistribution(mdf_mat_mat,'DisplayName','From ODF');
            %         hold off
            %         legend('-dynamicLegend','Location','northwest') % update legend
            %         saveFig(pfname, pname, name, [num2str(ii),'_gb_misor_dist_fit_', num2str(tolVal(ii)), 'deg']);
            %         % What we have plotted above is the uncorrelated misorientation angle
            %         % distribution for the Invar ODF. We can compare it to the uncorrelated
            %         % misorientation angle distribution of the uniform ODF
            
            %% Smooth of EBSD data
            if ~flagPlot
                F = splineFilter;
                ebsd_smoothed = smooth(ebsd(material),F,'fill',grains);
            end
            %         figure; set(gcf, 'Color', [1 1 1]);
            %         plot(ebsd_smoothed); hold on ;
            %         plot(grains.boundary,'linewidth',1.5); legend off;
            %         saveFig(pfname, pname, name, [num2str(ii),'_ebsd_smoothed_gb_', num2str(tolVal(ii)), 'deg']);
            %
            %% Orientations
            % compute mis2mean for the interpolated orientations
            if ~flagPlot
                [~,~,ebsd_smoothed.mis2mean] = calcGrains(ebsd_smoothed,'angle',10*degree);
            end
            
            % plot mis2mean for all phases
            %oM = ipdfHSVOrientationMapping(ebsd_smoothed(material).CS,ebsd_smoothed(material).CS);
            %         oM = TSLOrientationMapping(ebsd_smoothed(material).CS,ebsd_smoothed(material).CS);
            %         oM.maxAngle = tolVal(ii)*degree;
            %         color = oM.orientation2color(ebsd_smoothed(material).mis2mean);
            %         figure; set(gcf, 'Color', [1 1 1]);
            %         plot(ebsd_smoothed(material),color);
            %         hold on; plot(gB,'linewidth',1.5);
            %         legend off
            %         saveFig(pfname, pname, name, [num2str(ii), '_ebsd_smoothed_ori_', num2str(tolVal(ii)), 'deg']);
            %
            % Plot orientation map
            %oM = ipdfHSVOrientationMapping(ebsd(material));
            %clear oM;
            oM = TSLOrientationMapping(ebsd(material)); %clear color;
            color = oM.orientation2color(ebsd(material).orientations);
            figure; set(gcf, 'Color', [1 1 1]);
            plot(ebsd(material),color);
            saveFig(pfname, pname, name, [num2str(ii),'_ebsd_ori_', num2str(tolVal(ii)), 'deg']);
            
            % Plot orientation map
            clear oM;
            oM = TSLOrientationMapping(grains); clear color;
            color = oM.orientation2color(grains.meanOrientation);
            figure; set(gcf, 'Color', [1 1 1]);
            plot(grains, color);
            saveFig(pfname, pname, name, [num2str(ii),'_ebsd_oriMean_', num2str(tolVal(ii)), 'deg']);
            
            %         figure; set(gcf, 'Color', [1 1 1]);
            %         plot(grains, color);
            %         text(grains,int2str(grains.id));
            %         saveFig(pfname, pname, name, [num2str(ii),'_ebsd_oriMean_id_', num2str(tolVal(ii)), 'deg']);
            %
            % Plot IPF
            figure; set(gcf, 'Color', [1 1 1]);
            plot(oM);
            saveFig(pfname, pname, name, [num2str(ii),'_ori_ipf_', num2str(tolVal(ii)), 'deg']);
            
            %% Directions
            %         figure; set(gcf, 'Color', [1 1 1]);
            %         plot(grains,color,'micronbar','off');
            %
            %         % next we want to visualize the direction of the 100 axis
            %         dir = grains.meanOrientation * Miller(1,0,0,grains.CS);
            %
            %         % the lenght of the vectors should depend on the grain diameter
            %         len = 0.005*grains.diameter;
            %
            %         % arrows are plotted using the command quiver. We need to switch of auto
            %         % scaling of the arrow length
            %         hold on
            %         quiver(grains,len.*dir,'autoScale','off','color','black');
            %         hold off
            %         saveFig(pfname, pname, name, [num2str(ii),'_grains_dirArr_', num2str(tolVal(ii)), 'deg']);
            %
            %% Plot PDF / IPDF
            if ~flagPlot
                ori = ebsd_smoothed(material).orientations;
                odf = calcODF(ori,'Resolution', '5*degree');
                % try to compute an optimal kernel
                psi = calcKernel(ebsd_smoothed(material).orientations);
                psiMean = calcKernel(grains(material).meanOrientation);
                %compute the ODF with the kernel psi
                odf = calcODF(ori,'kernel',psi);
                %h = [Miller(1,0,0,odf.CS),Miller(1,1,0,odf.CS),Miller(1,1,1,odf.CS)];
                h = Miller(0,0,1,odf.CS);
            end
            
            %         figure; set(gcf, 'Color', [1 1 1]);
            %         plotPDF(odf,h,'colorrange','equal','antipodal','silent');
            %         mtexColorbar;
            %         mtexColorbar; % remove colorbars
            %         CLim(gcm,'equal');
            %         mtexColorbar; % add a single colorbar
            %         saveFig(pfname, pname, name, [num2str(ii),'_pdf_', num2str(tolVal(ii)), 'deg']);
            %
            figure; set(gcf, 'Color', [1 1 1]);
            plotIPDF(odf,h,'colorrange','equal','antipodal','silent');
            mtexColorbar;
            CLim(gcm, [0 3]);
            saveFig(pfname, pname, name, [num2str(ii),'_ipdf_', num2str(tolVal(ii)), 'deg']);
            %
            %         % define a unimodal ODF with the same modal orientation
            %         if ~flagPlot
            %             odf_model = unimodalODF(calcModes(odf),'halfwidth',15*degree);
            %         end
            %         figure; set(gcf, 'Color', [1 1 1]);
            %         plotPDF(odf_model,h,'antipodal','silent');
            %         mtexColorbar
            %         CLim(gcm,'equal');
            %         saveFig(pfname, pname, name, [num2str(ii),'_pdf_err_', num2str(tolVal(ii)), 'deg']);
            %         % compute the difference
            %         if ~flagPlot
            %             odfModel_odf = calcError(odf_model,odf);
            %         end
            %
            %         % % Zero Range Method
            %         if ~flagPlot
            %             odf_zero = calcODF(ori,'zero_range');
            %         end
            %         figure; set(gcf, 'Color', [1 1 1]);
            %         plotPDF(odf_zero,h,'antipodal','silent');
            %         mtexColorbar
            %         saveFig(pfname, pname, name, [num2str(ii),'_pdf_zeroRange_', num2str(tolVal(ii)), 'deg']);
            %
            %         if ~flagPlot
            %             pf = calcPoleFigure(odf,Miller(1,0,0,odf.CS),...
            %                 equispacedS2Grid('points',1000,'antipodal'));
            %         end
            %         figure; set(gcf, 'Color', [1 1 1]);
            %         %plot(pf,'MarkerSize',4);
            %         plot(pf,'contourf');
            %         mtexColorbar;
            %         %CLim(gcm,'equal');
            %         saveFig(pfname, pname, name, [num2str(ii),'_contourf_', num2str(tolVal(ii)), 'deg']);
            %
            %% Plot grain size, diameter and aspect ratios distributions
            % ids of the outer boundary segment
            if ~flagPlot
                outerBoundary_id = any(grains.boundary.grainId==0,2);
                % next we compute the corresponding grain_id
                grain_id = grains.boundary(outerBoundary_id).grainId;
                % remove all zeros
                grain_id(grain_id==0) = [];
                grainsCleaned = grains;
                grainsCleaned(grain_id) = [];
            end
            
            figure; set(gcf, 'Color', [1 1 1]);
            hist(grainsCleaned); legend off;
            grainsStat.meanArea = mean(grainsCleaned.area);
            grainsStat.maxArea = max(grainsCleaned.area);
            grainsStat.minArea = min(grainsCleaned.area);
            grainsStat.stdArea = std(grainsCleaned.area);
            saveFig(pfname, pname, name, [num2str(ii),'_grains_area_', num2str(tolVal(ii)), 'deg']);
            
            if ~flagPlot
                [d] = diameter(grainsCleaned);
            end
            figure; set(gcf, 'Color', [1 1 1]);
            hist(grainsCleaned.diameter); legend off;
            xlabel('Grains diameter (micron)');
            ylabel('Frequency (%)');
            grainsStat.meanDiam = mean(d);
            grainsStat.maxDiam = max(d);
            grainsStat.minDiam = min(d);
            grainsStat.stdDiam = std(d);
            saveFig(pfname, pname, name, [num2str(ii),'_grains_diam_', num2str(tolVal(ii)), 'deg']);
            
            % Aspect ratios
            figure; set(gcf, 'Color', [1 1 1]);
            % the size of the dots corresponds to the area of the grains
            scatter(grainsCleaned.shapeFactor, grainsCleaned.aspectRatio, 70*grainsCleaned.area./max(grains.area));
            xlabel('Length of the minor axis (micron)');
            ylabel('Length of the major axis (micron)');
            saveFig(pfname, pname, name, [num2str(ii),'_grains_shapeRatio_dist_', num2str(tolVal(ii)), 'deg']);
            
            figure; set(gcf, 'Color', [1 1 1]);
            plot(grainsCleaned,log(grainsCleaned.aspectRatio));
            mtexColorbar
            saveFig(pfname, pname, name, [num2str(ii),'_grains_shapeRatio_', num2str(tolVal(ii)), 'deg']);
            
            %% KAM;
            % ignore grain boundary misorientations
            if ~flagPlot
                kam_0 = ebsd_smoothed.KAM;
            end
            figure; set(gcf, 'Color', [1 1 1]);
            plot(ebsd_smoothed,kam_0 ./ degree);
            mtexColorbar
            saveFig(pfname, pname, name, [num2str(ii),'_KAM_', num2str(tolVal(ii)), 'deg']);
            
            % ignore misorientation angles > threshold
            if ~flagPlot
                kam_1 = KAM(ebsd_smoothed,'threshold',tolGB(ii)*degree);
            end
            figure; set(gcf, 'Color', [1 1 1]);
            plot(ebsd_smoothed,kam_1./degree);
            mtexColorbar
            saveFig(pfname, pname, name, [num2str(ii),'_KAM_thres_', num2str(tolVal(ii)), 'deg']);
            
            % consider also second order neigbors
            if ~flagPlot
                kam_2 = KAM(ebsd_smoothed,'order',2);
            end
            figure; set(gcf, 'Color', [1 1 1]);
            plot(ebsd_smoothed,kam_2./degree);
            mtexColorbar
            saveFig(pfname, pname, name, [num2str(ii),'_KAM_2ndOrder_', num2str(tolVal(ii)), 'deg']);
            
            % ignore misorientation angles > threshold && consider also second order neigbors
            if ~flagPlot
                kam_3 = KAM(ebsd_smoothed,'order',2,'threshold',tolGB(ii)*degree);
            end
            figure; set(gcf, 'Color', [1 1 1]);
            plot(ebsd_smoothed,kam_3./degree);
            mtexColorbar
            saveFig(pfname, pname, name, [num2str(ii),'_KAM_2ndOrder_thres_', num2str(tolVal(ii)), 'deg']);
            
            %% Save results
            res(ii).color = color;
            res(ii).CS = CS;
            res(ii).d = d;
            res(ii).ebsd = ebsd;
            res(ii).ebsd_smoothed = ebsd_smoothed;
            res(ii).ext = ext;
            res(ii).F = F;
            res(ii).filepath = filepath;
            res(ii).fname = fname;
            res(ii).gB = gB;
            res(ii).grain_id = grain_id;
            res(ii).grains = grains;
            res(ii).grainsCleaned = grainsCleaned;
            res(ii).grainsStat = grainsStat;
            res(ii).h = h;
            res(ii).id = id;
            res(ii).Indexed = Indexed;
            res(ii).kam_0 = kam_0;
            res(ii).kam_1 = kam_1;
            res(ii).kam_2 = kam_2;
            res(ii).kam_3 = kam_3;
            res(ii).mAngle = mAngle;
            res(ii).material = material;
            %        res(ii).mdf_mat_mat = mdf_mat_mat;
            res(ii).name = name;
            res(ii).notIndexed = notIndexed;
            res(ii).odf = odf;
            %         res(ii).odf_mat = odf_mat;
            %         res(ii).odf_model = odf_model;
            %         res(ii).odfModel_odf = odfModel_odf;
            %         res(ii).odf_zero = odf_zero;
            res(ii).ori = ori;
            res(ii).outerBoundary_id = outerBoundary_id;
            %         res(ii).pf = pf;
            res(ii).pfname = pfname;
            res(ii).pname = pname;
            res(ii).psi = psi;
            res(ii).psiMean = psiMean;
            res(ii).tol.tolVal = tolVal;
            res(ii).tol.tolGB = tolGB;
            res(ii).toRemove = toRemove;
            
        end
        save(fullfile([pname, '\results.mat']), 'res');
        display('EBSD data analysis finished !');
    else
        display('No EBSD results found in the given path !');
    end
end