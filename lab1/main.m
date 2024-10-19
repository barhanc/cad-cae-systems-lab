function main()
    img = imread('img.jpg');
    img = flip(img, 1);

    imshow(img)

    axis xy
    axis on
    grid on
    set(gca, 'xtick', 0:25:1150)

    hold on

    p = 3;
    knot_vector = cat(2, repelem([0], 1, p), 0:25:1150, repelem([1150], 1, p));
    coeffs = [325, 324, 324, 333, 350, 366, 373, 400, 415, 397, 374, 377, 313, 312, 305, 312, 324, 304, 278, 280, 277, 286, 289, 300, 296, 293, 295, 298, 285, 272, 265, 270, 272, 277, 279, 281, 284, 288, 288, 277, 290, 281, 280, 310, 310, 332, 341, 341, 356];
    precision = 100;
    splines_comp(precision, knot_vector, coeffs);
end
