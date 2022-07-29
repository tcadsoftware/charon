
SET (${tcad-charon_ENABLE_DAKOTA_DRIVERS} FALSE)

IF (${tcad-charon_ENABLE_DAKOTA_DRIVERS})

   execute_process (COMMAND echo RSYNCING Dakota Drivers)

   IF(IS_DIRECTORY src/Dakota)
       execute_process (COMMAND rm -fr src/Dakota)
   ENDIF()

   execute_process (COMMAND rsync -av  --exclude build  --progress ${CMAKE_SOURCE_DIR}/scripts/Dakota  src/ --delete)

   execute_process (COMMAND find ./src/Dakota/driver \( -name "*pyc" \) -exec rm {} \;)
   execute_process (WORKING_DIRECTORY src/Dakota/driver COMMAND python3  -m compileall ./)
   execute_process (WORKING_DIRECTORY src/Dakota/ COMMAND chmod -R g+rx driver)

   execute_process (WORKING_DIRECTORY src/ COMMAND ln -f -s Dakota/driver/driver.py charonDriver.py )
   execute_process (WORKING_DIRECTORY src/ COMMAND ln -f -s Dakota/driver/thresholdVoltage.py thresholdVoltage )

   #################################
   # install the Dakota Drivers upon make install
   #################################

   INSTALL(CODE "execute_process (WORKING_DIRECTORY src/Dakota COMMAND python3 installDakotaDrivers.cmake.py ${CMAKE_BINARY_DIR} ${CMAKE_INSTALL_PREFIX})")

   #################################
   # install the Dakota Drivers upon make install -- DONE
   #################################

ENDIF()


