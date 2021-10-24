# Unit Tests

This directory houses all of Charon's unit tests.  To add a new unit test,
simply use the `addUnitTest` script:
```
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ Usage:  addUnitTest TestName                                                ┃
┃                                                                             ┃
┃ This script will create a new group of unit tests named TestName and then   ┃
┃ add them into the CTest harness.  If you are creating the unit tests for a  ┃
┃ particular class in Charon, the recommendation would be to use the class    ┃
┃ name as TestName.                                                           ┃
┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
```

Note that there is a `UnitTestTemplate/` directory.  It's contents are used as
a template by the `addUnitTest` script to set up a new group of tests.  If at
any point we wish to modify what gets created by the `addUnitTest` script, we
need to modify the contents of that directory.
