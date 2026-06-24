function sessints = getSession_ints( tracking )
% return intervals for session
    sessints = struct();
    sessints.rest1 = [ tracking.events.subSessions(1,2)  tracking.events.subSessions(2,1) ];
    sessints.rest2 = [ tracking.events.subSessions(2,2)  tracking.events.subSessions(3,1) ];
    sessints.pre_mazeints = [ tracking.events.onMaze_startTime{1} tracking.events.openBox_ints{1}(:,2) ];
    sessints.pre_homeints = [ [ tracking.events.subSessions(1,1) ; tracking.events.openBox_ints{1}(1:end-1,2) ] tracking.events.openBox_ints{1}(:,1) ];
    sessints.learn_mazeints = [ tracking.events.onMaze_startTime{2} cellfun(@(x) x(end), tracking.events.rewardTime ) ];
    sessints.learn_homeints = [ [ tracking.events.subSessions(2,1) ; tracking.events.openBox_ints{2}(1:end-1,2) ] tracking.events.openBox_ints{2}(:,1) ];
    sessints.post_mazeints = [ tracking.events.onMaze_startTime{3} tracking.events.openBox_ints{3}(:,2) ];
    sessints.post_homeints = [ [ tracking.events.subSessions(3,1) ; tracking.events.openBox_ints{3}(1:end-1,2) ] tracking.events.openBox_ints{3}(:,1) ];

end