# HIBF

[![build status][1]][2]
[![codecov][3]][4]
[![license][5]][6]
![platforms][9]
<!-- [![latest release][7]][8] -->

<!--
    Above uses reference-style links with numbers.
    See also https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet#links.

    For example, `[![build status][1]][2]` evaluates to the following:
        `[link_text][2]`
        `[2]` is a reference to a link, i.e. `[link_text](https://...)`

        `[link_text]` = `[![build status][1]]`
        `[1]` is once again a reference to a link - this time an image, i.e. `[![build status](https://...)]
        `![build status]` is the text that should be displayed if the linked resource (`[1]`) is not available

    `[![build status][1]][2]` hence means:
    Show the picture linked under `[1]`. In case it cannot be displayed, show the text "build status" instead.
    The picture, or alternative text, should link to `[2]`.
-->

[1]: https://img.shields.io/github/actions/workflow/status/seqan/Hierarchical_Interleaved_Bloomfilter/ci_linux.yml?branch=main&style=flat&logo=github&label=CI "Open GitHub actions page"
[2]: https://github.com/seqan/Hierarchical_Interleaved_Bloomfilter/actions?query=branch%3Amain
[3]: https://codecov.io/gh/seqan/Hierarchical_Interleaved_Bloomfilter/branch/main/graph/badge.svg?token=BH1FQiBBle "Open Codecov page"
[4]: https://codecov.io/gh/seqan/Hierarchical_Interleaved_Bloomfilter
[5]: https://img.shields.io/badge/license-BSD-green.svg "Open Copyright page"
[6]: https://github.com/seqan/Hierarchical_Interleaved_Bloomfilter/blob/main/LICENSE.md
[7]: https://img.shields.io/github/release/seqan/Hierarchical_Interleaved_Bloomfilter.svg "Get the latest release"
[8]: https://github.com/seqan/Hierarchical_Interleaved_Bloomfilter/releases/latest
[9]: https://img.shields.io/badge/platform-linux%20%7C%20bsd%20%7C%20osx-informational.svg

This library contains the HIBF and layout algorithm.

## Quick start

To use the HIBF lib in your app:

```cmake
include (FetchContent)
FetchContent_Declare (
    hibf_fetch_content
    GIT_REPOSITORY "https://github.com/seqan/Hierarchical_Interleaved_Bloomfilter"
    GIT_TAG "main")
option (INSTALL_HIBF "" OFF)
FetchContent_MakeAvailable (hibf_fetch_content)

# ...

target_link_libraries (<your_app> PUBLIC seqan::hibf)
```

## Sponsorships

[![Vercel](https://raw.githubusercontent.com/seqan/Hierarchical_Interleaved_Bloomfilter/main/test/documentation/.vercel/powered-by-vercel.svg)](https://vercel.com/?utm_source=seqan&utm_campaign=oss)

Vercel is kind enough to sponsor our documentation preview-builds within our pull requests. Check them out!
