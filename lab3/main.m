function main()

filename="apple.jpg";
elementsx=4;
elementsy=4;
maxerror=1;
color_edges_black=true;
mode="b";
 
n_frames = 8;
for iter=1:n_frames
    disp(iter);
    [im, totel] = bitmap_h(filename, elementsx, elementsy, maxerror, iter, color_edges_black, mode);
    
    fig = figure;
    imshow(im);
    title(sprintf("Iter: %d Total elements: %d", iter, totel));
    fig = getframe(fig);
    
    imwrite(fig.cdata, sprintf('%s_frame_%d.jpg',mode,iter));
end
end