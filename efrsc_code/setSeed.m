%% Sets random number generator seed
function setSeed(seed)
  if (isOctave)
    rand('state',seed);
  else 
    if exist('RandStream')
      RandStream.setDefaultStream(RandStream('mt19937ar','seed',seed));
    else
      rand('twister',seed);
    end
  end
end
