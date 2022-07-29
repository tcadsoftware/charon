IF(NOT IS_DIRECTORY src/interpreter)
       execute_process (COMMAND mkdir src/interpreter)
ENDIF()
IF(IS_DIRECTORY src/interpreter/parsers)
        execute_process (COMMAND rm -fr src/interpreter/parsers)
        execute_process (COMMAND rm -fr src/interpreter/modifiers)
ENDIF()
execute_process (COMMAND rm -fr src/interpreter )
execute_process (COMMAND rsync -av  --exclude generateParserDocumentation --exclude parsers --exclude modifiers --progress ${CMAKE_SOURCE_DIR}/scripts/charonInterpreter  src/interpreter --delete)

#execute_process (WORKING_DIRECTORY 

execute_process (WORKING_DIRECTORY src/interpreter/charonInterpreter/parseGenerator COMMAND python3 generateInterpreter.py --verbosity=20)
execute_process (COMMAND find ./src/interpreter \( -name "*pyc" \) -exec rm {} \;)
execute_process (COMMAND python3 --version)
execute_process (COMMAND mv src/interpreter/charonInterpreter/parseGenerator parseGenerator)
execute_process (WORKING_DIRECTORY src/interpreter/charonInterpreter/ COMMAND python3 -m compileall ./)
execute_process (WORKING_DIRECTORY src/ COMMAND chmod -R g+rX ./)
execute_process (WORKING_DIRECTORY src/ COMMAND ln -f -s interpreter/charonInterpreter/charonInterpreter.py charonInterpreter )
execute_process (WORKING_DIRECTORY src/ COMMAND ln -f -s interpreter/charonInterpreter/charonInterpreter.py charonInterpreter.py )
execute_process (WORKING_DIRECTORY src/ COMMAND ln -f -s interpreter/charonInterpreter/charonInterpreter.py chirp )
execute_process (WORKING_DIRECTORY src/ COMMAND ln -f -s interpreter/charonInterpreter/tools/compareParameterLists.py compareParameterLists )
execute_process (WORKING_DIRECTORY src/ COMMAND ln -f -s interpreter/charonInterpreter/tools/createInputTest.py createInputTest )
execute_process (WORKING_DIRECTORY src/ COMMAND ln -f -s interpreter/charonInterpreter/tools/dagify.py dagify )
execute_process (WORKING_DIRECTORY src/ COMMAND ln -f -s interpreter/charonInterpreter/tools/xmlToLCM.py xmlToLCM )
execute_process (COMMAND mv parseGenerator src/interpreter/charonInterpreter/)

#################################
# install the interpeter upon make install
#################################
  INSTALL(CODE "execute_process (WORKING_DIRECTORY src/interpreter/charonInterpreter COMMAND python3 installInterpreter.cmake.py ${CMAKE_BINARY_DIR} ${CMAKE_INSTALL_PREFIX})")

#################################
# install the interpeter upon make install -- DONE
#################################

