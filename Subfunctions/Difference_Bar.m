function indMaxBand = Difference_Bar(Signal_d1,Signal_d2,fAxis,ylimit,yShift,color)

thrNeighborPts = 1;
thrFreq = 0;
Signal_1 = Signal_d1;
Signal_2 = Signal_d2;
Diff = (Signal_1 - Signal_2)>= 0;
freqAxis = fAxis;
cc = bwconncomp(Diff);
if(cc.NumObjects~=0)
    numPixels = cellfun(@numel,cc.PixelIdxList);
    cc.PixelIdxList = cc.PixelIdxList(numPixels>=thrNeighborPts);
    % Threshold lower bound for frequency
    cc.PixelIdxList = cc.PixelIdxList(freqAxis(cellfun(@min,cc.PixelIdxList))>=thrFreq);
    % Cluster with the highest area under
    sig_all = cellfun(@(x) Signal_1(x),cc.PixelIdxList,'UniformOutput',0);
    thr_all = cellfun(@(x) Signal_2(x),cc.PixelIdxList,'UniformOutput',0);
    [~,indMaxCluster] = max(cellfun(@sum,cellfun(@(x,y) x-y,sig_all,thr_all,'UniformOutput',0)));
    %                                 [~,indMaxCluster] = (cellfun(@sum,cellfun(@(x,y) x-y,sig_all,thr_all,'UniformOutput',0)));
    
    if(~isempty(indMaxCluster))
        for i =1:size(cc.PixelIdxList,2)
            indMaxBand = cc.PixelIdxList{i};
            %             if (all(indMaxBand>=15&indMaxBand<=30))
            %                 GridChans_sgnf_15_30(i) =1;
            %             else
            %                 GridChans_sgnf_15_30(i) =0;
            %             end
            freqBand = freqAxis(indMaxBand);
            if size(freqBand,2)>=2
            AreaUnder_PLV_Prc = sum(Signal_1(indMaxBand)-Signal_2(indMaxBand));
            semilogx(freqAxis(indMaxBand),yShift*ones(length(freqAxis(indMaxBand)),1)','Color',color,'LineWidth',2)
            end
            %             ylim([ylimit(1),ylimit(2)]);
                significant_bar{i} = freqAxis(indMaxBand);
                hold on;
            
        end
    else
        indMaxBand = 0;
        freqBand = [];
        AreaUnder_PLV_Prc = 0;
        AreaUnder_PLV_Prc_SelectedBand = 0;
%         significant_bar = [];
    end
    indMaxBand =   significant_bar;
    %     indMaxBand = GridChans_sgnf_15_30;
else
    indMaxBand = 0;
end