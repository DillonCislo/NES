function compile_minimize_elastic_energy

mex -v -O minimize_elastic_energy.cpp -I/usr/include:/usr/local/include -L/usr/lib:/usr/local/lib -lCGAL -lgmp -lboost_thread
end

