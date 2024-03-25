using DelimitedFiles

function readIpoptOutput(folderName, fileName)
    content = readdlm(folderName*fileName)
    iters = content[end-16,4]
    time = content[end-1,6]
    return iters, time
end