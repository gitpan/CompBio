# MySQL dump 8.14
#
# Host: localhost    Database: PDB
#--------------------------------------------------------
# Server version	3.23.41

#
# Table structure for table 'blast'
#

DROP TABLE IF EXISTS Blast;
CREATE TABLE Blast (
  ID INT UNSIGNED ZEROFILL NOT NULL,
  Match_ID INT UNSIGNED ZEROFILL NOT NULL,
  EffectiveIdentities FLOAT(8,3) UNSIGNED NOT NULL,
  EquivalentIdentities TINYINT UNSIGNED NOT NULL,
  Qstart SMALLINT UNSIGNED DEFAULT NULL,
  Qend SMALLINT UNSIGNED DEFAULT NULL,
  Mstart SMALLINT UNSIGNED DEFAULT NULL,
  Mend SMALLINT UNSIGNED DEFAULT NULL,
  PRIMARY KEY (ID,Match_ID)
) TYPE=MyISAM PACK_KEYS=1 CHECKSUM=1 COMMENT="Corrilation using Blast";

#
# Table structure for table 'Annotation'
#

DROP TABLE IF EXISTS Annotation;
CREATE TABLE Annotation (
  ID INT UNSIGNED ZEROFILL NOT NULL,
  Definition TEXT NOT NULL,
  Supplemental TEXT,
  Submitted TEXT,
  PRIMARY KEY  (ID)
) TYPE=MyISAM PACK_KEYS=1 CHECKSUM=1 COMMENT="Original and extra annotation";


#
# Table structure for table 'keywords'
#

DROP TABLE IF EXISTS Keywords;
CREATE TABLE Keywords (
  ID INT UNSIGNED ZEROFILL NOT NULL,
  GeneName VARCHAR(12) DEFAULT NULL,
  Function VARCHAR(32) DEFAULT NULL,
  PRIMARY KEY  (ID),
  KEY (GeneName)
) TYPE=MyISAM PACK_KEYS=1 CHECKSUM=1 COMMENT="Single term fields, such as HosH or amidotransferase";

#
# Table structure for table 'Peptide'
#

DROP TABLE IF EXISTS Peptide;
CREATE TABLE Peptide (
  ID INT UNSIGNED ZEROFILL NOT NULL,
  Sequence TEXT NOT NULL,
  PRIMARY KEY  (ID)
) TYPE=MyISAM PACK_KEYS=1 CHECKSUM=1;
