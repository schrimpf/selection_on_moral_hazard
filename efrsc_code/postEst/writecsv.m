function writecsv(filename,colnames,data)
  out = fopen(filename,'w');
  for i=1:numel(colnames)
    fprintf(out,'%s,',colnames{i});
  end
  fprintf(out,'\n');
  for i=1:size(data,1)
    fprintf(out,'%.16g,',data(i,:));
    fprintf(out,'\n');
  end
  fclose(out);
end