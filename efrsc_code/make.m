%% Compiles c-code -- has very limited make-like functionality
% (modification date checking). Set force to true if you want to rebuild
% regardless of modification date.
function make(force)
  % we need to link to nlopt and mpfr. adjust these variables so
  % that those libraries can be found. Many linux distributions include mpfr,
  % so you might not need to change anything for it. If you install nlopt
  % with default options, it will be /usr/local and this will work. 
  local = '/usr/local';
  Lib = ['-L' local '/lib'];     % add to linktime library path
  Inc = ['-I' local '/include']; % add into include path
  linkLib = ['''-Wl,-rpath,' local '/lib''']; % add to runtime library path 

  % options needed to link to lapack with standard fortran interface 
  % you will need to modify this to fit your system and preferred
  % version of lapack. You can use the version of lapack included with matlab
  % or octave. Here we use AMD's math core library, but it makes little
  % difference. For octave, it is likely that -llapack suffices, for matlab
  % -lmwlapack likely would work. Note that with these libraries, you should
  % change chol.c and sampleMu.c to use the fortran interface to lapack (see
  % comments in those files)
  lapack = ['-L/opt/acml5.1.0/gfortran64_fma4_mp/lib -lacml_mp ' ...
            '''-Wl,-rpath=/opt/acml5.1.0/gfortran64_fma4_mp/lib'' ' ...
            '-lgfortran ' ... 
            '-I/opt/acml5.1.0/gfortran64_fma4_mp/include'];
  lapack = '-llapack';
  if (~isOctave)  % Matlab gets its compilation options from mexopts.sh, not
                  % environmental variables directly, so write a version of
                  % mexopts.sh that respects the putenv() commands below. You
                  % may need to modify this for your system.
    warning('Overwriting mexopts.sh');
    file = fopen('mexopts.sh','w');
    fprintf(file,['    TMW_ROOT="$MATLAB"\n' ...
                  'MFLAGS=''''\n' ...
                  'if [ "$ENTRYPOINT" = "mexLibrary" ]; then \n' ...
                  'MLIBS="-L$TMW_ROOT/bin/$Arch -lmx -lmex -lmat -lmwservices ' ...
                  '-lut -lm"\n'  ...
                  'else  \n'...
                  'MLIBS="-L$TMW_ROOT/bin/$Arch -lmx -lmex -lmat -lm" \n' ...
                  'fi \n\n' ...
                  'CC="$CC"\n' ...
                  'CXX="$CXX"\n' ...
                  'FC="$F77"\n' ...
                  'RPATH="-Wl,-rpath-link,$TMW_ROOT/bin/$Arch"\n' ... 
                  'CFLAGS="$CFLAGS -fno-omit-frame-pointer -pthread"\n' ...
                  'CLIBS="$RPATH $MLIBS -lm -lstdc++"\n' ...
                  'CXXFLAGS="$CXXFLAGS -fPIC -fno-omit-frame-pointer -pthread"\n' ...
                  'CXXLIBS="$CLIBS"\n' ...
                  'LD="$CC" \n' ...
                  'LDEXTENSION=''.mexglx''\n' ... % change to mexa64 on 64 bit system
                  'LDFLAGS="-pthread -shared -Wl,--version-script,$TMW_ROOT/' ...
                  'extern/lib/$Arch/$MAPFILE -Wl,--no-undefined"\n']); 
    fclose(file);    
    putenv = @(name,value) setenv(name,value);
  end

  %% set compilers and options, change if needed 
  putenv('CC','gcc');
  putenv('CXX','g++');
  putenv('F77','gfortran');
  putenv('CFLAGS',['-O0 -fopenmp -std=c99 -fPIC ' ...
		   '-march=native -Wall -g']); % change march as needed for
                                             % your system. If you get
                                             % SIGILL: illegal
                                             % instruction, you have the
                                             % wrong march. Using 'uname -m'
                                             % the output will work, but you
                                             % might get better results by
                                             % being more specific. native
                                             % will work for newer versions
                                             % of gcc
  putenv('CPPFLAGS','');
  putenv('CXXFLAGS',['-D_GLIBCXX_DLL -march=native -O2 -fPIC ' ...
		     '-Wall -ffast-math -fopenmp ']);
  llink = ['-lnlopt -lmpfr -lgmp -lpthread -lgomp'];
  if (nargin<1)
    force = false;
  end
  if (needCompile('alcoa','.o') || force)
    eval(sprintf('mex -v -c -lmpfr -lgmp %s %s alcoa.c', ...
                 Inc,Lib));
    checkCompile('alcoa','.o');
  end
  
  if (~exist('arms','dir')) % try to download Adaptive Rejection Metropolis
                            % Sampling (ARMS)
    [dl txt] = system(['wget http://www1.maths.leeds.ac.uk/~wally.gilks/' ...
                       'adaptive.rejection/arms.method/arms_method.zip']) 
    fprintf('%s\n',txt);
    if (dl~=0) 
      fprintf(['ERROR: failed to download ARMS. Perhaps the link is ' ...
               'dead.\n']);
      exit(dl);
    end
    [zip txt]=system('unzip arms_method.zip -d arms');
  end
  
      
  if (needCompile('arms/arms','.o') || force);
    mex -v -c arms/arms.c -o arms/arms.o;     
    if (~isOctave) % matlab ignores the -o option
      system('mv arms.o arms/.');
    end
    checkCompile('arms/arms','.o');
  end
  if(isOctave)
    ext = '.mex';
  else 
    arch = getenv('ARCH');
    if (arch == 'glnx86')
      ext = '.mexglx'; 
    else
      ext = '.mexa64'; % assumes 64 bit linux
    end
  end 
  if (needCompile('XtTimesKronSIdTimesY',ext) || force)
    mex XtTimesKronSIdTimesY.c -v;
    checkCompile( 'XtTimesKronSIdTimesY',ext);
  end
    
  if (isOctave)     
    % octave's builtin chol hangs unexpectedly. rebuild it from lapack.
    if (needCompile('chol','_c.mex') || force)
      eval(sprintf('mex -v %s chol.c -o chol_c ',lapack));
      checkCompile('chol','_c.mex');
    end
  end
  if(isOctave)
    ext = '_c.mex';
  else 
    arch = getenv('ARCH');
    if (arch == 'glnx86')
      ext = '_c.mexglx'; 
    else
      ext = '_c.mexa64'; % assumes 64 bit linux
    end
  end 
  for f={'sampleLambdaOmega','sampleLambdaOmegaMultMH','samplePsi','findChoices', ...
         'findChoicesNMH', 'sampleSigL','findValidLatent3', ...
         'randdtn','sampleShape', 'sampleLamlo', ...
         'findChoicesMultMH', 'findChoicesNmultMH'} 
    file = f{1};
    if (needCompile(file,ext) || force) 
      cmd=sprintf(['mex  -v '...
                   ' -Iarms arms/arms.o alcoa.o ' ...
                   '%s %s %s  %s' ...
                   ' %s.c -o %s_c'],Inc,Lib,llink, linkLib,file,file);
      fprintf('%s\n',cmd);
      eval(cmd);
      checkCompile(file,ext);
    end
  end
  file = 'sampleMu';       
  if (needCompile(file,ext) || force)
    % this is separate because needs linking with lapack
    cmd=sprintf(['mex  -v '...
                 ' -Iarms arms/arms.o alcoa.o ' ...
                 '%s %s %s  %s %s ' ...
                 ' %s.c -o %s_c'],Inc,Lib,llink, linkLib,lapack,file,file);
    fprintf('%s\n',cmd);
    eval(cmd);
    checkCompile(file,ext);
  end
    
  
end

% check if we need to (re) compile  file.c into fileext
function n=needCompile(file,ext) 
  dm = dir('make.m');
  dh = dir([file '.c']);
  d2 = dir([file ext]);
  n=(isempty(d2) || dh.datenum>d2.datenum || ...
        dm.datenum>d2.datenum);
end

% check if we successfully compiled file.c into fileext
function checkCompile(file,ext) 
  d1 = dir([file '.c']);
  check = dir([file ext]);
  if (isempty(check) || check.datenum<d1.datenum)
    fprintf('##########################################\n');
    error(['Failed to compile ' file]);	      
  else 
    fprintf(['\nSuccessfully compiled %s\n' ...
             '##################################\n\n' ...
            ],file);
  end
end

% function addToMexOpts(name,value)
%   file = fopen('mexopts.sh','a');
%   fprintf(file,'%s=''%s'' \n',name,value);
%   fclose(file);  
% end
