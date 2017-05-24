function Build-Unix
{
    bash.exe -c 'mkdir -p build-unix && cd build-unix && cmake --target all --config Debug -- -j 10 .. && make'
}

function Run-Unix
{
    Build-Unix
    bash.exe -c 'echo && ./bin/swsh && echo'
}