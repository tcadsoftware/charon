TRIBITS_PROJECT_DEFINE_EXTRA_REPOSITORIES(
  Trilinos        Trilinos             GIT git@cee-gitlab.sandia.gov:Charon/Trilinos        "HASPACKAGES,POST" Continuous
  src             src                  GIT git@cee-gitlab.sandia.gov:Charon/src             "HASPACKAGES,POST" Continuous
  nightlyTests    test/nightlyTests    GIT git@cee-gitlab.sandia.gov:Charon/nightlyTests    "NOPACKAGES"       Continuous
  nightlyTestsOUO test/nightlyTestsOUO GIT git@cee-gitlab.sandia.gov:Charon/nightlyTestsOUO "NOPACKAGES"       Continuous
  heavyTests      test/heavyTests      GIT git@cee-gitlab.sandia.gov:Charon/heavyTests      "NOPACKAGES"       Continuous
  heavyTestsOUO   test/heavyTestsOUO   GIT git@cee-gitlab.sandia.gov:Charon/heavyTestsOUO   "NOPACKAGES"       Continuous
  docs            docs                 GIT git@cee-gitlab.sandia.gov:Charon/docs            "NOPACKAGES"       Continuous
  )
