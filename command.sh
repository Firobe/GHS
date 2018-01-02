./gengraph cycle 37 -format dot | dot -Tdot | sed 's/^\s//g' | sed 's/ //g' | sed 's/\s/ /g' > pouet.dot 
