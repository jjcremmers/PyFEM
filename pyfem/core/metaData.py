# SPDX-License-Identifier: MIT
# Copyright (c) 2011-2026 Joris J.C. Remmers

"""Analysis metadata helpers."""

from datetime import datetime, timezone
from os import environ
from subprocess import DEVNULL, CalledProcessError, check_output

from pyfem import __version__
from pyfem.util.dataStructures import Properties
from pyfem.util.logger import getLogger, separator


def get_github_build() -> str:
    """Return GitHub Actions build metadata or the local Git commit.

    :returns: GitHub build identifier, local Git commit hash, or ``"unknown"``.
    """
    github_run = environ.get("GITHUB_RUN_NUMBER")
    github_attempt = environ.get("GITHUB_RUN_ATTEMPT")
    github_sha = environ.get("GITHUB_SHA")

    if github_run and github_sha:
        build = f"github-run-{github_run}"
        if github_attempt:
            build += f".{github_attempt}"
        return f"{build} ({github_sha[:12]})"

    try:
        return check_output(
            ["git", "rev-parse", "HEAD"],
            stderr=DEVNULL,
            text=True,
        ).strip()
    except (CalledProcessError, FileNotFoundError):
        return "unknown"


def get_wall_clock_time() -> datetime:
    """Return the current local wall-clock time.

    :returns: Timezone-aware local wall-clock timestamp.
    """
    return datetime.now(timezone.utc).astimezone()


def store_analysis_metadata(globdat, wall_clock_time: datetime | None = None) -> None:
    """Store wall-clock and build metadata on the global data object.

    :param globdat: Global data object receiving the ``metadata`` attribute.
    :param wall_clock_time: Wall-clock timestamp captured at CLI startup.
    """
    if wall_clock_time is None:
        wall_clock_time = get_wall_clock_time()

    globdat.metadata = Properties(
        {
            "wallClockTime": wall_clock_time.isoformat(timespec="seconds"),
            "githubBuild": get_github_build(),
            "pyfemVersion": __version__,
        }
    )


def print_analysis_metadata(globdat) -> None:
    """Print analysis metadata stored on the global data object.

    :param globdat: Global data object containing a ``metadata`` attribute.
    """
    metadata = globdat.metadata
    logger = getLogger()

    separator("-")
    logger.info("  Analysis metadata")
    logger.info("    Wall clock time ......... : %s", metadata.wallClockTime)
    logger.info("    GitHub build ............ : %s", metadata.githubBuild)
    logger.info("    PyFEM version ........... : %s", metadata.pyfemVersion)
    separator("-")
