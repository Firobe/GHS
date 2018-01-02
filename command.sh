./gengraph planar 75 3 3 -format dot | dot -Tdot | sed 's/^\s//g' | sed 's/ //g' | sed 's/\s/ /g' > pouet.dot 

