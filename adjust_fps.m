function adjust_fps(v_in,v_out,fps_out)

vid = VideoReader(v_in);
writerObj = VideoWriter(v_out,'mpeg-4');
writerObj.FrameRate = fps_out;  % Change the desired frame rate here.
open(writerObj);
nFrames = vid.NumberOfFrames;
vidHeight = vid.Height;
vidWidth = vid.Width;
mov(1:nFrames) = struct('cdata', zeros(vidHeight, vidWidth, 3, 'uint8'),'colormap', []);
for k = 1:length(mov)
    mov(k).cdata = read(vid, k);
    writeVideo(writerObj, mov(k).cdata);
end
% hf = figure;
% set(hf, 'position', [20 150 vidWidth vidHeight]);
% movie(hf, mov, 1, vid.FrameRate);
close(writerObj);