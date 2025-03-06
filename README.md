# nbbp

## Overview

`nbbp` is an R package for inferring the effective reproduction number and dispersion parameter from stuttering chains.

## Getting started

The R package can be installed with `remotes::install_github("CDCGov/nbbp")`, to install just the package, or with `remotes::install_github("CDCGov/nbbp", build_vignettes = TRUE)` to install the package and build vignettes.

(You can also clone or download this repository, `cd` to the resulting folder, and install the package with `devtools::install()` or, to build the vignettes, `devtools::install(build_vignettes = TRUE)`.)

Some of the package's dependencies require additional system-level dependencies to work, such as [pandoc](https://pandoc.org/).
The [Posit Package Manager](https://packagemanager.posit.co/client/#/) is useful for determining what needs to be run to install these on your operating system, e.g. [for rstan](https://packagemanager.posit.co/client/#/repos/cran/packages/overview?search=rstan).

When installing, the package will compile a [stan model](https://mc-stan.org/users/interfaces/rstan) and print a fair bit of text to screen.
This is expected and can be ignored as long as the installation is successful.

If building the vignettes, after loading the package in R (`library("nbbp")`), the following vignettes will be available:

- `vignette("nbbp")`: introduction to the package and inference when all chains are completely observed and extinct
- `vignette("advanced_data")`: how to handle censored observations and very large (or non-extinct) chains
- `vignette("default_priors")`: plots of default priors and prior predictive distributions

## Project Admin

Andy Magee, PhD, (@afmagee42)

---

---

## General Disclaimer

This repository was created for use by CDC programs to collaborate on public health related projects in support of the [CDC mission](https://www.cdc.gov/about/organization/mission.htm). GitHub is not hosted by the CDC, but is a third party website used by CDC and its partners to share information and collaborate on software. CDC use of GitHub does not imply an endorsement of any one particular service, product, or enterprise.

## Public Domain Standard Notice

This repository constitutes a work of the United States Government and is not
subject to domestic copyright protection under 17 USC § 105. This repository is in
the public domain within the United States, and copyright and related rights in
the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).
All contributions to this repository will be released under the CC0 dedication. By
submitting a pull request you are agreeing to comply with this waiver of
copyright interest.

## License Standard Notice

This repository is licensed under ASL v2 or later.

This source code in this repository is free: you can redistribute it and/or modify it under
the terms of the Apache Software License version 2, or (at your option) any
later version.

This source code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the Apache Software License for more details.

You should have received a copy of the Apache Software License along with this
program. If not, see http://www.apache.org/licenses/LICENSE-2.0.html

The source code forked from other open source projects will inherit its license.

## Privacy Standard Notice

This repository contains only non-sensitive, publicly available data and
information. All material and community participation is covered by the
[Disclaimer](https://github.com/CDCgov/template/blob/master/DISCLAIMER.md)
and [Code of Conduct](https://github.com/CDCgov/template/blob/master/code-of-conduct.md).
For more information about CDC's privacy policy, please visit [http://www.cdc.gov/other/privacy.html](https://www.cdc.gov/other/privacy.html).

## Contributing Standard Notice

Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo)
and submitting a pull request. (If you are new to GitHub, you might start with a
[basic tutorial](https://help.github.com/articles/set-up-git).) By contributing
to this project, you grant a world-wide, royalty-free, perpetual, irrevocable,
non-exclusive, transferable license to all users under the terms of the
[Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or
later.

All comments, messages, pull requests, and other submissions received through
CDC including this GitHub page may be subject to applicable federal law, including but not limited to the Federal Records Act, and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

## Records Management Standard Notice

This repository is not a source of government records but is a copy to increase
collaboration and collaborative potential. All government records will be
published through the [CDC web site](http://www.cdc.gov).

## Additional Standard Notices

Please refer to [CDC's Template Repository](https://github.com/CDCgov/template) for more information about [contributing to this repository](https://github.com/CDCgov/template/blob/main/CONTRIBUTING.md), [public domain notices and disclaimers](https://github.com/CDCgov/template/blob/main/DISCLAIMER.md), and [code of conduct](https://github.com/CDCgov/template/blob/main/code-of-conduct.md).
