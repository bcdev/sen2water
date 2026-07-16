#
# Copyright (C) 2025 ESA
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import pytest
from click.testing import CliRunner

from eopf import __version__
from eopf.cli.cli import eopf_cli

from eopf.config.config import EOConfiguration


# @pytest.mark.parametrize("opts", [("--help",)])
# @pytest.mark.parametrize("opts", [
#     #("trigger","local","/home/martin/projects/sen2water/eopf/tests/test_payload.yaml")
#      ("trigger","local","/home/hannes/projects/opt-mpc/sen2water/eopf/tests/test_resamplingpu_payload.yaml")
# ])
@pytest.mark.parametrize("opts", [("trigger","local","E:/olaf/bc/workspace/sen2water_reengineering-experiment-resampling-pu/sen2water/eopf/tests/test_idepixclassificationpu_payload.yaml")])
# @pytest.mark.parametrize("opts", [("trigger","local","E:/olaf/bc/workspace/sen2water_reengineering-experiment-resampling-pu/sen2water/eopf/tests/test_resampling_payload.yaml")])
# @pytest.mark.parametrize("opts", [("trigger","local","E:/olaf/bc/workspace/sen2water_reengineering-experiment-resampling-pu/sen2water/eopf/tests/test_payload.yaml")])
# @pytest.mark.parametrize("opts", [("trigger","local",".\tests\test_payload.yaml")])
def test_cli(opts):
    # set env var EOPF_CONFIGURATION_FOLDER to absolute path of directory that contains file eopf.toml
    EOConfiguration().load_default_file()
    runner = CliRunner()
    r = runner.invoke(eopf_cli, args=opts, catch_exceptions=False)
    print(r.stderr)
    print(r.output)
    assert r.exit_code == 0
    # fmt: off

