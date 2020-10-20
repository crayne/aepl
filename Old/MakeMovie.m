%%  SegIms  -  Epl
%       part of Epl
%       little to no changes made
%

function MakeMovie(Fr,Name)

[P,F,~] = fileparts(Name);

P2 = fullfile(P,'MP4s');
if ~exist(P2,'dir')
    mkdir(P2)
end 

VidName = fullfile(P2,[F,'.mp4']);

v = VideoWriter(VidName,'MPEG-4');

v.FrameRate = 15;  % Default 30
v.Quality = 50;    % Default 75
open(v);

for t = 1:size(Fr,2)
 writeVideo(v,Fr{t});
end

close(v);
end