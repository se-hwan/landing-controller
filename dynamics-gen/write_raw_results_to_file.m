function write_raw_results_to_file(data, filename)
 
fid = fopen(filename,'w');
fwrite(fid,data,'single');

end