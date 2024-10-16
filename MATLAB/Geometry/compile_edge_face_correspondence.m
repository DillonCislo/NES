function compile_edge_face_correspondence

clc;

[projectDir, ~, ~] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(projectDir);

CGALDir = '/mnt/data0/Code/Cpp/CGAL-6.0/include';
boostDir = '/usr/include/boost';
eigenDir = '/mnt/data0/Code/Cpp/eigen-3.4.0';
iglDir = '/mnt/data0/Code/Cpp/libigl/include';
tbbDir = '/usr/include/tbb';

CXXOPTIMFLAGS = '"-O3" ';
CXXFLAGS = '"$CXXFLAGS -march=native -fopenmp -fPIC -std=c++17" ';
LDFLAGS = '"$LDFLAGS -fopenmp" ';
CGALCompileFlags = '-DCGAL_EIGEN3_ENABLED=true ';

includeFlags = [ ... 
    '-I' eigenDir ' ' ... 
    '-I' iglDir ' ' ... 
    '-I' boostDir ' ' ... 
    '-I' CGALDir ' ' ... 
    '-I' tbbDir ' ' ... 
    '-I/usr/include:/usr/local/include ' ];

libFlags = '-L/usr/lib:/usr/local/lib:/usr/lib/x86-64-linux-gnu ';
libFlags = [ libFlags '-lgmp -lmpfr -lboost_thread -lboost_system -ltbb' ];

mexString = [ 'mex -R2018a -v -O edge_face_correspondence.cpp ' ... 
    'CXXOPTIMFLAGS=' CXXOPTIMFLAGS ' ' ... 
    'CXXFLAGS=' CXXFLAGS ' ' ... 
    'LDFLAGS=' LDFLAGS ' ' ... 
    CGALCompileFlags includeFlags libFlags ];

eval(mexString);

end

