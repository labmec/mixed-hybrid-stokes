file(GLOB headers *.h json.hpp) # creates a set headers with all the .h files
file(GLOB sources *.cpp) # creates a set sources with all the .cpp files
set(stokessources ${headers} ${sources} PARENT_SCOPE)