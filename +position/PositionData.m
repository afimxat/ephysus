classdef PositionData < neuro.basic.ChannelTimeData & ...
        matlab.mixin.indexing.RedefinesParen
    %POSITIONDATA Class for representing positional data in a 3D space over time.
    % This class stores the X, Y, Z coordinates and time information, and provides
    % methods for manipulating and accessing the position data.

    properties
        units % Units of measurement for the positional data (e.g., 'cm', 'meters').
        Info % Additional information about the position data.
    end

    methods
        function obj = PositionData(X,Y,Z,time)
            %PositionData Construct an instance of this class
            % Inputs:
            %   X, Y, Z: Arrays of coordinates in their respective dimensions.
            %   time: A TimeIntervalCombined object representing the time points for the data.
            % The constructor initializes a PositionData object either by copying another PositionData object's properties,
            % or by creating a new object with specified X, Y, Z, and time data, validating the input sizes.
            if nargin>0
                if isa(X, 'optiTrack.PositionData')
                    obj.data = X.data;
                    obj.time = X.timeIntervalCombined;
                    obj.units = X.units;
                    return;
                end

                if numel(X) ~= numel(Y) || numel(Z) ~= numel(Y) || numel(X) ~= time.getNumberOfPoints()
                    error('Sizes of XYZ or time are not equal.');
                end

                data = [X(:), Y(:), Z(:)];
                obj.data = array2table(data, 'VariableNames', {'X', 'Y', 'Z'});
                obj.time = time;
                obj.units = 'cm';
            end
        end

        function [positionData, idx] = getWindow(obj, range)
            % GETWINDOW Returns a subset of the PositionData within a specified time range.
            % Inputs:
            %   range: A 2-element vector [startTime endTime] specifying the time range.
            % Outputs:
            %   positionData: A new PositionData object for the specified time range.
            %   idx: Indices of the data points within the specified time range.
            ticd = obj.time;
            positionData = obj;
            window = ticd.getTimeIntervalForTimes(range);
            positionData.time = window;
            samples = ticd.getSampleForClosest(range);
            idx = samples(1):samples(2);
            positionData.data = obj.data(idx,:);
        end
        function obj=plus(obj,pd)
            % Overloads the plus operator to concatenate the data from another PositionData object
            % to the current object, ensuring they are sequential in time and have matching sample rates.
            if obj.time.getEndTime>=pd.time.getStartTime
                error(['Position data starts(%s) before the ' ...
                    'original ends(%s).'],pd.time.getStartTime, ...
                    obj.time.getEndTime)
            end
            if obj.time.getSampleRate~=pd.time.getSampleRate
                error('Samplerates are different. First:%dHz, Second:%dHz', ...
                    obj.time.getSampleRate,pd.time.getSampleRate)
            end
            newtime=obj.time+pd.time;
            newdata=[obj.data; pd.data];
            obj.time=newtime;
            obj.data=newdata;
        end
        function data=getData(obj)
            % Returns the data table of the PositionData object.
            data=obj.data;
        end

        function obj=setData(obj,data)
            % Sets the data table of the PositionData object.
            obj.data=data;
        end

        function mat=flatten2(obj)
            % Returns the first two rows of the data table, effectively "flattening" the data into 2D.
            mat=obj.data(1:2,:);
        end

        function mat=flatten3(obj)
            % Returns the first three rows of the data table, keeping the data in 3D.
            mat=obj.data(1:3,:);
        end
        function pdman = getManifold(obj)
            % Creates a manifold representation of the positional data and saves it to a file.
            % This method applies dimensionality reduction and manifold learning techniques on the data.
            % Refer to the method body for detailed operations including plotting and file saving.

            time = obj.time;
            timestr = matlab.lang.makeValidName(time.tostring);
            file1 = java.io.File(obj.source);
            manifoldFile = fullfile(char(file1.getParent), ...
                ['position.PositionDataManifold' timestr '.mat']);

            if exist(manifoldFile, 'file')
                s = load(manifoldFile);
                pdman = s.pdman;
                if isempty(pdman.units)
                    pdman.units = 'cm';
                end
                return;
            else
                manifoldFile1 = fullfile(char(file1.getParent), ...
                    ['position.PositionDataManifold' '*' '.mat']);
                fs=dir(manifoldFile1);
                if numel(fs)==1
                    s = load(fullfile(fs.folder,fs.name));
                    pdman = s.pdman;
                    if isempty(pdman.units)
                        pdman.units = 'cm';
                    end
                    return;
                end
            end

            try
                close(123);
            catch
            end
            figure(123);
            f = gcf;
            f.Position(3:4) = [2500 1500];
            tiledlayout(2, 3);

            c.numberOfPoints = 300;
            c.neighbors = 7;

            sampleRate = obj.time.getSampleRate;
            obj1 = obj.getDownsampled(sampleRate * 0.1);
            spd = obj1.getSpeed(3).Values;
            [~, I] = sort(spd, 'descend', 'MissingPlacement', 'last');

            data1 = table2array(obj1.data)';
            data2 = data1(:, I(1:(size(data1, 2) / 20)));
            data2(:, any(isnan(data2))) = [];

            nexttile(1, [2 1]);
            obj.plot3DtimeContinuous;
            title('Original');

            nexttile(3, [1 1]);
            manifold = external.Manifold.Manifold("Description");
            manifold = manifold.createGraph(data2, 'verbose', 'neighbors', c.neighbors, 'numPoints', c.numberOfPoints);
            manifold.plotGraph;
            ax = gca;
            ax.DataAspectRatio = [1 1 1];
            title('Graph');

            manifold = manifold.shortestPath('verbose');

            nexttile(6, [1 1]);
            manifold = manifold.scale('plot', 'sammon');
            ax = gca;
            ax.DataAspectRatio = [1 1 1];
            title('Scaled');

            pdman = position.PositionDataManifold(obj, manifold);
            pdman.config = c;

            nexttile(2, [2 1]);
            pdman.plot3DtimeContinuous;
            title('Dimension Reduced');

            ff = logistics.FigureFactory.instance(char(file1.getParent));
            ff.save(['position-PositionDataManifold-' timestr]);
            save(manifoldFile, 'pdman', '-mat');
        end
        function [velocity] = getSpeed(obj, smoothingWindowInSeconds)
            % Computes the speed from the positional data.
            % Inputs:
            %   smoothingWindowInSeconds: The window size for smoothing speed data.
            % Outputs:
            %   velocity: A Channel object representing the speed at each time point.
            data = table2array(obj.getData)';
            timeDiffSeconds = diff(seconds(obj.time.getTimePoints));
            timeDiffSeconds2=[timeDiffSeconds median(timeDiffSeconds)];
            squaredDiffs = zeros(size(data));
            for dimIndex = 1:size(data, 1)
                squaredDiffs(dimIndex, 1:(end-1) )= diff(data(dimIndex, :)).^2;
            end
            speeds = sqrt(sum(squaredDiffs, 1))./timeDiffSeconds2;
            if exist('smoothingWindowInSeconds', 'var')
                speeds = smoothdata(speeds, 'gaussian', obj.time.getSampleRate * ...
                    smoothingWindowInSeconds);
            end
            velocity = neuro.basic.Channel('Velocity', speeds, obj.time);
        end

        function [om]= getOccupancyMap(obj,xedges,zedges)
            % Generates an occupancy map from the positional data.
            % Inputs:
            %   xedges, zedges: Bin edges for the X and Z dimensions.
            % Outputs:
            %   om: An OccupancyMap object representing the spatial occupancy.
            if nargin==1
                om=neuro.placeField.OccupancyMap(obj,obj.time.getSampleRate);
            else
                om=neuro.placeField.OccupancyMap(obj, ...
                    obj.time.getSampleRate,xedges,zedges);
            end
            om.Units=obj.units;
        end
        function ax = plot(obj)
            % Plots the positional data over time.
            % This method generates a plot of the position data, automatically downsampling if necessary.
            numPointsInPlot = 100000;
            time = obj.time;
            t_org = seconds(time.getTimePoints() - (time.getZeitgeberTime() - time.getStartTime()));
            downsampleFactor = max(1, round(numel(t_org) / numPointsInPlot)); % add lower bound
            data = table2array(obj.getData());
            downsampledData = downsample(data, downsampleFactor); % remove loop
            t = hours(seconds(downsample(t_org, downsampleFactor)));
            plot(t, downsampledData);
            legend(obj.getData().Properties.VariableNames);
            xlabel('ZT (Hrs)');
            ylabel(['Location (', obj.units, ')']);
            ax = gca;
        end

        function ax = plot2D(obj, numPointsInPlot)
            % Generates a 2D scatter plot of the position data, downsampling if specified.
            % Input:
            %   numPointsInPlot: The number of points to include in the plot, for downsampling.
            % Set default value for numPointsInPlot if not provided
            if nargin < 2
                numPointsInPlot = 10000;
            end

            % Get time data
            ticd = obj.time;
            t_org = ticd.getTimePointsZT;

            % Calculate downsample factor
            downsampleFactor = round(numel(t_org) / numPointsInPlot);

            % Get input data
            dims = {'X','Z'};
            inputData = table2array(obj.data(:, dims));

            % Downsample input data
            downsampledData = downsample(inputData, downsampleFactor)';

            % Get colormap
            colorMap = linspecer(size(downsampledData, 2));

            % Plot scatter plot
            scatter(downsampledData(1,:), downsampledData(2,:), [], colorMap, ...
                'filled', 'MarkerFaceAlpha', .2, 'MarkerEdgeAlpha', .2, 'SizeData', 5);

            % Set axis labels and aspect ratio
            ylabel([dims{2} ' ' obj.units]);
            xlabel([dims{1} ' ' obj.units]);
            ax = gca;
            ax.DataAspectRatio = [1 1 1];

            % Set colormap
            colormap(colorMap);
        end
        function p = plot2DContinuous(obj, numPointsInPlot)
            % Generates a 2D plot of the position data with continuous lines, downsampling if specified.
            % Input:
            %   numPointsInPlot: The number of points to include in the plot, for downsampling.
            if ~exist('numPointsInPlot','var')
                numPointsInPlot=10000;
            end

            % Get data and time information
            data = obj.getData();
            time = obj.time;
            t_org = time.getTimePointsZT();

            % Downsample data
            downsampleFactor = round(numel(t_org) / numPointsInPlot);
            if downsampleFactor < 1
                downsampleFactor = 1;
            end
            data = downsample(data, downsampleFactor);
            t = downsample(t_org, downsampleFactor);

            % Plot continuous lines
            dims = {'X','Z'};
            p = plot(data.(dims{1}), data.(dims{2}));
            p.LineWidth = 2;
            p.Color = 'k';
            ylabel([dims{2} ' ' obj.units]);
            xlabel([dims{1} ' ' obj.units]);
        end
        function ax = plot3Dtime(obj, numPointsInPlot)
            % Generates a 3D scatter plot of the position data over time, downsampling if specified.
            % Input:
            %   numPointsInPlot: The number of points to include in the plot, for downsampling.
            % Set default value for numPointsInPlot if not provided
            if nargin < 2
                numPointsInPlot = 10000;
            end

            % Get time vector and downsample factor
            ticd = obj.time;
            t_org = seconds(ticd.getTimePointsZT);
            downsampleFactor = max(round(numel(t_org) / numPointsInPlot), 1);

            % Get data and downsample along the first two dimensions
            dims = {'X', 'Z', 'Y'};
            data0 = obj.getData();
            data1 = table2array(data0(:, dims))';
            data2 = zeros(size(data1));
            for ich = 1:size(data1, 1)
                data2(ich, :) = downsample(medfilt1(data1(ich,:), ...
                    ticd.getSampleRate()), downsampleFactor);
            end

            % Add time to the third dimension
            data2(3, :) = data2(3, :) + linspace(1, t_org(end) - t_org(2), size(data2, 2));

            % Plot data as a 3D scatter plot
            color1 = linspecer(size(data2, 2));
            scatter3(data2(1,:), data2(2,:), data2(3,:), [], color1, ...
                'filled', 'MarkerFaceAlpha', .2, 'MarkerEdgeAlpha', .2, 'SizeData', 5);
            xlabel(dims{1});
            ylabel(dims{2});
            zlabel('Time (s)');
            ax = gca;
            ax.DataAspectRatio = [1 1 4];
        end

        function ax = plot3DtimeContinuous(obj, numPointsInPlot)
            % Generates a 3D plot of the position data over time with continuous lines, downsampling if specified.
            % Input:
            %   numPointsInPlot: The number of points to include in the plot, for downsampling.
            if ~exist('numPointsInPlot','var')
                numPointsInPlot=10000;
            end
            colors=colororder;
            ticd=obj.time;
            t_org=seconds(ticd.getTimePointsZT);
            downsamplefactor=round(numel(t_org)/numPointsInPlot);
            if downsamplefactor<1
                downsamplefactor=1;
            end

            dims={'X','Z','Y'};
            data0=obj.getData;
            data1=table2array(data0(:,dims))';
            for ich=1:size(data1,1)
                data2(ich,:)=downsample(data1(ich,:),downsamplefactor); %#ok<AGROW>
            end
            lintime=linspace(0,t_org(end)-t_org(2),size(data2,2));
            data2(3,:)=data2(3,:)+lintime; % add time
            plot3(data2(1,:),data2(2,:),data2(3,:),LineWidth=1,Color= ...
                colors(1,:));
            hold on
            zlabel('Time (s)');
            ylabel(dims{2});
            xlabel(dims{1});
            ax=gca;
            ax.DataAspectRatio=[1 1 4];
            ax.ZDir="reverse";
            ax.YDir="reverse";
            %             nanidx=any(isnan(data2));
            %             nanarr=nan(size(data2));
            %             nanarr(1:2,nanidx)=0;
            %             nanarr(3,nanidx)=lintime(nanidx);
            %             plot3(nanarr(1,:),nanarr(2,:),nanarr(3,:),
            % LineWidth=2,Color=colors(2,:));
            %             data2(1:2,:)=0;
            %             data2(3,~nanidx)=lintime(~nanidx);
            %             plot3(data2(1,:),data2(2,:),data2(3,:),
            % LineWidth=2,Color=colors(3,:));
        end
        function ax = plot3D(obj, numPointsInPlot)
            % Generates a 3D scatter plot of the position data, downsampling if specified.
            % Input:
            %   numPointsInPlot: The number of points to include in the plot, for downsampling.
            if ~exist('numPointsInPlot','var')
                numPointsInPlot=10000;
            end
            ticd=obj.time;
            t_org=ticd.getTimePointsZT;
            downsamplefactor=round(numel(t_org)/numPointsInPlot);

            dims={'X','Z','Y'};
            data0=obj.getData;
            data1=table2array(data0(:,dims))';
            for ich=1:size(data1,1)
                data2(ich,:)=downsample(medfilt1(data1(ich,:), ...
                    ticd.getSampleRate),downsamplefactor); %#ok<AGROW>
            end
            plot3(data2(1,:),data2(2,:),data2(3,:));
            zlabel(dims{3});
            ylabel(dims{2});
            xlabel(dims{1});
        end
        function [] = plot3DMark(obj,mark,color)
            % Plots markers on a 3D plot of the position data at specified indices.
            % Inputs:
            %   mark: Indices of the data points to mark.
            %   color: Color of the markers.
            if ~exist('color','var')
                color=[];
            end
            ticd=obj.time;
            t_org=seconds(ticd.getTimePoints);

            dims={'X','Z','Y'};
            data0=obj.getData;
            data1=table2array(data0(:,dims))';
            data1(3,:)=data1(3,:)+linspace(1,t_org(end)-t_org(2), ...
                size(data1,2)); % add time

            data2=data1(:,mark);
            try
                s=scatter3(data2(1,:),data2(2,:),data2(3,:),[],color,'filled', ...
                    'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5,'SizeData',30);
            catch ME
                s=scatter3(data2(1,:),data2(2,:),data2(3,:),'filled', ...
                    'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5,'SizeData',30);
            end
            if isempty(color)
                s.MarkerFaceColor="k";
            end
            ylabel([dims{2} ' ' obj.units]);
            xlabel([dims{1} ' ' obj.units]);

        end
        function [] = plot2DMark(obj,mark)
            % Plots markers on a 2D scatter plot of the position data at specified indices.
            % Input:
            %   mark: Indices of the data points to mark.
            dims={'X','Z','Y'};
            data0=obj.getData;
            data1=table2array(data0(:,dims))';

            data2=data1(:,mark);
            s=scatter(data2(1,:),data2(2,:),[],'filled', ...
                'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5,'SizeData',5);
            s.MarkerFaceColor="k";
            ylabel([dims{2} ' ' obj.units]);
            xlabel([dims{1} ' ' obj.units]);

        end
        function obj = getTimeWindow(obj,timeWindow)
            % Selects a subset of the PositionData within a specified time window.
            % Input:
            %   timeWindow: A 2-element vector [startTime endTime] specifying the time window.
            ticd=obj.time;
            ticdnew=ticd.getTimeIntervalForTimes(timeWindow);
            s1=ticd.getSampleForClosest(ticdnew.getStartTime);
            s2=ticd.getSampleForClosest(ticdnew.getEndTime);
            obj.time=ticdnew;
            obj.data=obj.data(s1:s2,:);
        end
        function obj = getDownsampled(obj,dsfactor)
            % Returns a downsampled version of the PositionData object.
            % Input:
            %   dsfactor: Downsample factor.
            data1=table2array(obj.data)';
            for idim=1:size(data1,1)
                data1ds(idim,:)=downsample(medfilt1(data1(idim,:),dsfactor), ...
                    dsfactor,dsfactor-1);
            end
            obj.data=array2table(data1ds',"VariableNames", ...
                obj.data.Properties.VariableNames);
            ticd=obj.time;
            obj.time=ticd.getDownsampled(dsfactor);
        end
        function obj = getMedianFiltered(obj,winseconds)
            % Applies a median filter to the positional data.
            % Input:
            %   winseconds: The window size for the median filter in seconds.
            data1=table2array(obj.data)';
            data2=smoothdata(data1,2,'movmedian',winseconds*obj.time.getSampleRate);
            obj.data=array2table(data2',"VariableNames", ...
                obj.data.Properties.VariableNames);
        end
        function obj = getMeanFiltered(obj,winseconds)
            % Applies a mean filter to the positional data.
            % Input:
            %   winseconds: The window size for the mean filter in seconds.
            data1=table2array(obj.data)';
            data2=smoothdata(data1,2,'movmean',winseconds*obj.time.getSampleRate);
            obj.data=array2table(data2',"VariableNames", ...
                obj.data.Properties.VariableNames);
        end
        function [data1,sample] = getPositionForTimes(obj,times)
            % Retrieves the position data for specified time points.
            % Input:
            %   times: An array of time points for which to retrieve data.
            % Outputs:
            %   data1: The position data for the specified time points.
            %   sample: The indices of the data points corresponding to the specified times.
            if numel(times)>0
                time=obj.time;
                sample=time.getSampleForClosest(times);
                data1=obj.data(sample,:);
            else
                sample=[];
                data1=obj.data(sample,:);
            end
        end

        function [obj, folder]= saveInPlainFormat(obj,folder,ext1)
            % Saves the PositionData object in a plain text format.
            % Inputs:
            %   folder: The directory to save the files in.
            %   ext1: The file extension for the position points data.
            if ~exist('ext1','var')
                ext1='position.points.csv';
            end
            extt='position.time.csv';
            if exist('folder','var')
                if ~isfolder(folder)
                    folder= pwd;
                    warning(['The folder does not exist. ' ...
                        'Files will be created in %s.'],folder);
                end
            else
                folder= fileparts(obj.source);
            end
            time=obj.time; %#ok<*PROPLC>
            timestr=matlab.lang.makeValidName(time.tostring);
            time.saveTable(fullfile(folder,[timestr extt]));
            file1=fullfile(folder,[timestr ext1]);
            writetable(obj.data,file1);
            file2=fullfile(folder,[timestr '.mat']);
            save(file2,'obj');
            folder=string(py.os.path.realpath(py.os.path.expanduser(folder)));
            %             obj=obj.loadPlainFormat(folder);
        end
        function obj= loadPlainFormat(obj,folder)
            % Loads the PositionData object from plain text format files.
            % Input:
            %   folder: The directory containing the files.
            ext1='position.points.csv';
            extt='position.time.csv';
            [file1, uni]=obj.getFile(folder,ext1);
            obj.source=file1;
            obj.data=readtable(obj.source);
            folder=fileparts(file1);
            obj.time=time.TimeIntervalCombined( ...
                fullfile(folder,[uni extt]));
        end
        function [file2, uni]=getFile(~,folder,extension)
            % Helper method to retrieve files based on a specified extension.
            % Inputs:
            %   folder: The directory to search in.
            %   extension: The file extension to search for.
            % Outputs:
            %   file2: The full path to the file found.
            %   uni: A unique identifier derived from the file name.
            if ~exist('folder','var')
                folder= pwd;
            end
            if isfile(folder)
                [folder1,name,ext1]=fileparts(folder);
                uni1=split([name ext1],extension);
                uni=uni1{1};
                file1=dir(fullfile(folder1,[uni,extension]));
            else
                file1=dir(fullfile(folder,['*' extension]));
                if numel(file1)>1
                    [name,folder1] = uigetfile({['*' extension],extension}, ...
                        'Selectone of the position files',folder);
                    file1=dir(fullfile(folder1,name));
                end
            end
            file2=fullfile(file1.folder,file1.name);
            [~,name,ext1]=fileparts(file2);
            uni1=split([name ext1],extension);
            uni=uni1{1};
        end
    end
    methods (Access=protected)
        function obj = parenReference(obj, indexOp)
            % Overrides the default behavior for referencing elements using parentheses.
            try
                idx=indexOp.Indices{:};
            catch ME
                idx=indexOp(1,1).Indices{:};
            end
            data=table2array(obj.data);
            data(~idx,:) =nan;
            obj.data=array2table(data,VariableNames= ...
                obj.data.Properties.VariableNames);
        end

        function obj = parenAssign(obj,indexOp,varargin)
            % Overrides the default behavior for assigning elements using parentheses.
            if isempty(obj)
                obj = varargin{1};
            end
            if isscalar(indexOp)
                assert(nargin==3);
                rhs = varargin{1};
                obj.ContainedArray.(indexOp) = rhs.ContainedArray;
                return;
            end
            [obj.(indexOp(2:end))] = varargin{:};
        end

        function n = parenListLength(obj,indexOp,ctx)
            % Determines the number of elements in a list during a parenthetical reference operation.
            if numel(indexOp) <= 2
                n = 1;
                return;
            end
            containedObj = obj.(indexOp(1:2));
            n = listLength(containedObj,indexOp(3:end),ctx);
        end

        function obj = parenDelete(obj,indexOp)
            % Overrides the default behavior for deleting elements using parentheses.
            obj.ContainedArray.(indexOp) = [];
        end
    end
    methods
        function out = cat(dim,varargin)
            % Overloads the cat function for concatenating arrays of this class.
            numCatArrays = nargin-1;
            newArgs = cell(numCatArrays,1);
            for ix = 1:numCatArrays
                if isa(varargin{ix},'ArrayWithLabel')
                    newArgs{ix} = varargin{ix}.ContainedArray;
                else
                    newArgs{ix} = varargin{ix};
                end
            end
            out = ArrayWithLabel(cat(dim,newArgs{:}));
        end

        function varargout = size(obj,varargin)
            % Overloads the size function to provide the size of the data array.
            [varargout{1:nargout}] = [1 1];
        end
    end
end

