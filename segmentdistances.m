clear all

a = sqrt(3)/4;
b = 1/4;
c = 1/2;
d = sqrt(3)/2;

last_layer = 8;
ncells = 3*last_layer*(last_layer+1)

%Partitioning the line segment
rball = 1.;
nh = 100 %number of segments horizontally
x_sep = rball/nh;
x_pts = 0:x_sep:rball;
n_boxes = nh;

%Vertices of hexagonal cells (U_i)
e = [a b; 0 0.5; -a b; -a -b; 0 -0.5; a -b];
%Vectors between centers of cells
dirs = [1 0; c d; -c d; -1 0; -c -d; c -d; 1 0];
%Points for distances from segments to nearest 6 cells
Points = [1 0; d c; 0 1; -d c; 0 -1; d -c];

%Compute nearest points to origin for each cell
dists = zeros(n_boxes,ncells);
layer = 2:last_layer;
tic
if true
  toc
  %Choose box
  for j = 1:n_boxes
    j
    xl = x_pts(j);
    xr = x_pts(j+1);

    %Set distance from first 6 cells as dist from segment to point in Points
    for i = 1:6
      Point = Points(i,:);

      cvx_begin quiet
        variables x(1,2);
        minimize norm( Point - x);
        subject to
          xl <= x(1) <= xr;
          x(2) == 0;
      cvx_end

      dists(j,i) = norm(Point - x);
    end

    %For further cells, use CVX to determine shortest distance between segment and cell
    cell_ind = 7;
    %Choose cell layer
    for l = layer
      rad = l*sqrt(3)/2; %Set distance of first cell in each layer/sector
      %Choose sector
      for sect = 1:6
        %Choose cell in sector, l cells per sector at layer l
        for i = 0:(l-1)
          %Cells in layer l and sector i will be between first cells in sectors i and i+1
          center = rad*(dirs(sect,:)*(1-i/l) + dirs(sect+1,:)*(i/l));

          %Minimize distance between x and hexagonal cell, t is coefficient to move in direction of cell's vertices e
          cvx_begin quiet
            variables t(1,6) x(1,2);
            minimize norm( t*e + center - x);
            subject to
              t >= 0;
              sum(t) == 1;
              norm(x) <= rball; %redundant but needed for earlier version
              xl <= x(1) <= xr;
              x(2) == 0;
          cvx_end

          dists(j,cell_ind) = norm(t*e + center - x);
          %Increment cell index to keep count
          cell_ind = cell_ind+1;
        end
      end
    end
    toc
  end

  %Save to file with date
  DT = datetime();
  Day = day(DT);
  Month = month(DT);
  Year = year(DT);
  filename = sprintf('partition_dists_%04d%02d%02d_npartitions_%d.mat', Year, Month, Day, nh);
  save(filename,'dists');
end
