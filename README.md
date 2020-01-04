# cellular-automata-for-PCP
http://49.233.23.19/dokuwiki/doku.php?id=projects:rotationprojects:developmentcodingofeleganse_yuanzhe_201911


For the project deacribtion please read the wiki page above.

CA.py --- Rule 1.0
CA_hex_none_divide.py ---Rule 1.2
CA_hex_none_divide_2.0.py ---Rule 2.0

----
Note for CA.py

The code can be simply divide to two parts, **class Playground** and **class Cell**.
  
For each Cell object, it got its **position(x,y)** and **other attributes**(TFs contains e.g.). And it got its basic function **divid()**, once know its neighbor, combined with the cell attributes, this function return two daughtercells(Cell type).
  
The Playground class handles the upgrading of each step by **hulahula()**, and to provide information for any action involving more than one cell. **update_state()** contains rules of CA.
  
For now, the only rule is if there is an empty space, then the cell divide. If not, it remains to the next upgrade.
----

CA_hex_none_divide.py
