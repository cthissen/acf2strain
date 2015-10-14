function dualfprintf(fid,message)
    fprintf(fid,message);
    fprintf(message);
end