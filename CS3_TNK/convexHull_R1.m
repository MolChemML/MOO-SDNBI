function [facetIndices, facetNormals] = convexHull_R1(z,dim)
    % Returns the facet indicies defining the facets, and the facet
    % normals. The last element of the facet normals is the offset (it isnt
    % part of the vector).
    % I am using the qhull C++ code as it returns the facet normals
    % (ratherr
    % than MATLAB). 
    %
    % The code creates an input 'data.txt', with format:
    %       first line contains the dimension
    %       second line contains the number of input points
    %       remaining lines contain point coordinates
    %
    % Then runs the function (qconvex) in the OS terminal that outputs an 'output.txt'
    % The structure of output text is defined by the optional arguments to
    % qconvex: i (print vertices incident to each facet)
    %       The first line is the number of facets. The remaining lines list the vertices for each facet, one facet per line. 
    %       The indices are 0-relative indices of the corresponding input points. 
    %       The facets are oriented. Option 'Fv' displays an unoriented list of vertices with a vertex count per line. Options 'o' and 'Ft' displays coordinates for each vertex prior to the vertices for each facet.
    %       Simplicial facets (e.g., triangles in 3-d) consist of d vertices. Non-simplicial facets in 3-d consist of 4 or more vertices. For example, a facet of a cube consists of 4 vertices. Use option 'Qt' to triangulate non-simplicial facets.
    %       For 4-d and higher convex hulls and 3-d and higher Delaunay triangulations, d vertices are listed for all facets. A non-simplicial facet is triangulated with its centrum and each ridge. The index of the centrum is higher than any input point. Use option 'Fv' to list the vertices of non-simplicial facets as is. Use option 'Ft' to print the coordinates of the centrums as well as those of the input points.
    %
    % Second input n (prints hyperplane normals with offsets)
    %       The first line is the dimension plus one. 
    %       The second line is the number of facets. 
    %       The remaining lines are the normals for each facet, one normal per line. The facet's offset follows its normal coefficients.
    %       The normals point outward, i.e., the convex hull satisfies Ax <= -b where A is the matrix of coefficients and b is the vector of offsets.
    %       A point is inside or below a hyperplane if its distance to the hyperplane is negative. A point is outside or above a hyperplane if its distance to the hyperplane is positive. Otherwise a point is on or coplanar to the hyperplane.
    %       If cdd output is specified ('FD'), Qhull prints the command line, the keyword "begin", the number of facets, the dimension (plus one), the keyword "real", and the normals for each facet. The facet's negative offset precedes its normal coefficients (i.e., if the origin is an interior point, the offset is positive). Qhull ends the output with the keyword "end".
    %
    % Then reads the output and converts to a matrix where each row
    % corresponds to a different facet.
  
    fclose('all');
    
    fout = sprintf('convHullOut2_%s.txt',['_', datestr(now,'mm-dd-yyyy-HH-MM-SS')]);
    fin  = sprintf('convHullIn2_%s.txt',['_',datestr(now,'mm-dd-yyyy-HH-MM-SS')]);
    
    fclose(fopen(fout, 'w')); % empty file for output 3
    
    % Format text file with data
    fname = fin;
    fileID = fopen(fname,'w');
    numberOfPoints  = size(z,1);
    string          = repmat('%6e ', 1, dim);
    string          = [string, '\r\n'];
    fprintf(fileID,'%d\n',dim);
    fprintf(fileID,'%d\n',numberOfPoints);
    fprintf(fileID, string, z');    
    fclose(fileID);
    
    % Run qconvex function

    system(join(['qconvex Fv n <', fin, '>', fout]));
     system(join(['qconvex GD3 <', fin, '>', 'eg.01.cube']));
    % Read output
    fid = fopen(fout);
    tline = fgetl(fid);
    numberOfFacets = str2double(tline);
    linecount = 1;
    facetIndices = cell(numberOfFacets,1);
    
    for i = 1:numberOfFacets
        tline = fgetl(fid);
        facetIndices{i} = str2num(tline(2:end)) + 1;
    end
    
    tline = fgetl(fid); % skip 2 lines (no. facets, no. indices)
    tline = fgetl(fid);
    
    facetNormals = [];
    for i = 1:numberOfFacets
        tline = fgetl(fid);
        facetNormals(end+1,:) = -str2num(tline); %facing inside convex hull
    end

    facetNormals = facetNormals(:,1:end-1);
    facetNormals = facetNormals;
    fclose(fid);
    fclose('all');
    
   % voronoi(
    
    delete(fin)
    delete(fout)
    
    
    
  