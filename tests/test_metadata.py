from zeizify.czi_metadata import (
    read_bit_ceiling,
    read_channel_names,
    read_field_names,
    read_pixel_size_um,
    read_scenes,
)

XML = """<?xml version="1.0"?>
<ImageDocument><Metadata>
  <Experiment><ExperimentBlocks><AcquisitionBlock>
    <!-- decoy: acquisition-track channel order differs from the image C axis -->
    <Channel Id="Channel:9" Name="DAPI" /><Channel Id="Channel:8" Name="AF594" /><Channel Id="Channel:7" Name="AF488" />
  </AcquisitionBlock></ExperimentBlocks></Experiment>
  <Information><Image>
    <ComponentBitCount>14</ComponentBitCount>
    <Dimensions>
      <Channels>
        <Channel Id="Channel:0" Name="AF488"><ComponentBitCount>14</ComponentBitCount></Channel>
        <Channel Id="Channel:1" Name="AF594" />
        <Channel Id="Channel:2" Name="DAPI" />
      </Channels>
      <S><Scenes>
        <Scene Index="0" Name="P1"><ArrayName>A1</ArrayName><Shape Name="A1"><ColumnIndex>0</ColumnIndex><RowIndex>0</RowIndex></Shape></Scene>
        <Scene Name="P2" Index="1"><ArrayName>A01</ArrayName><Shape Name="A1"><ColumnIndex>0</ColumnIndex><RowIndex>0</RowIndex></Shape></Scene>
        <Scene Index="2" Name="P1"><ArrayName>B12</ArrayName><Shape Name="B12"><ColumnIndex>11</ColumnIndex><RowIndex>1</RowIndex></Shape></Scene>
      </Scenes></S>
    </Dimensions>
  </Image></Information>
  <Scaling><Items>
    <Distance Id="X"><Value>3.25e-07</Value></Distance>
    <Distance Id="Y"><Value>3.25e-07</Value></Distance>
  </Items></Scaling>
</Metadata></ImageDocument>"""


def test_channels_come_from_information_block_not_first_seen():
    assert read_channel_names(XML) == ["AF488", "AF594", "DAPI"]


def test_scenes_parse_regardless_of_attribute_order():
    scenes = read_scenes(XML)
    assert [(s.index, s.name, s.array_name) for s in scenes] == [(0, "P1", "A1"), (1, "P2", "A01"), (2, "P1", "B12")]
    assert scenes[2].row_index == 1 and scenes[2].col_index == 11
    assert read_field_names(XML) == {0: "P1", 1: "P2", 2: "P1"}


def test_scaling_and_bits():
    assert abs(read_pixel_size_um(XML) - 0.325) < 1e-9
    assert read_bit_ceiling(XML) == 16383


def test_missing_values_are_none():
    xml = "<ImageDocument><Metadata><Information><Image/></Information></Metadata></ImageDocument>"
    assert read_pixel_size_um(xml) is None
    assert read_bit_ceiling(xml) is None
    assert read_scenes(xml) == []
    assert read_channel_names(xml) == []
