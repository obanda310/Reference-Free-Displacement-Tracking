classdef Trajectory < Experiment
    properties
        posIni   % position of trajectory in first frame
        posNow   % position of trajectory in current frame
        dPos     % posNow - posIni
        frameIni   % first frame this trajectory appears in
        frameEnd   % last frame this trajectory appears in
        
        rightNeighbor % trajectory to right of this trajectory
        
        row      % row that trajectory is a member of
        modeRowX   % mode of all round_dx0 trajectory values in this row > threshold
        modeRowY   % mode of all round_dx0 trajectory values in this row > threshold
        newX    % error-adjusted x value
        newY    % error-adjusted y value
        newdx0  % error-adjusted dx0 value
        newdy0  % error-adjusted dx0 value
    end
    methods
        function obj = Trajectory(frameNo,trajNo)
            if frameNo > obj.numFrames
                error('The frameNo is greater than the number of frames in this experiment!')
            elseif trajNo > obj.numTrajs
                error('The trajNo is greater than the number of trajectories in this experiment!')
            end
            
            thisInfo = obj.trajData((obj.trajData(:,2)==trajNo),:);
            obj.posIni = [thisInfo(1,4),thisInfo(1,5)];
            obj.posNow = [thisInfo(frameNo,4),thisInfo(frameNo,5)];
            obj.dPos = obj.posNow - obj.posIni;
            obj.frameIni = thisInfo(1,3)+1;
            obj.frameEnd = thisInfo(end,3)+1;
        end
    end
end