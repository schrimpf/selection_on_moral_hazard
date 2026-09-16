%% return true if in octave, false if matlab
function s = isOctave() 
  s = exist('OCTAVE_VERSION','builtin');
end