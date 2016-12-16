<?xml version="1.0"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
                version="1.0">

  <!-- Copyright 2014 The MathWorks, Inc. -->

  <xsl:output method="html" indent="no"/>

  <xsl:param name="demosroot"/>
  <xsl:param name="helpdir"/>
  <xsl:param name="exampledir"/>

  <xsl:template match="toc">
    <ul>
      <xsl:apply-templates select="tocitem"/>
    </ul>
    <xsl:if test="$exampledir != ''">
      <ul>
        <li id="demos.xml">
          <div class="node_element">
            <div class="wrapper img_wrapper"> </div>
            <div class="wrapper text_wrapper">
              <a href="{concat($demosroot,'/demos.xml')}" class="corrected_url">Examples</a>
            </div>
          </div>
        </li>
      </ul>
    </xsl:if>
  </xsl:template>

  <xsl:template match="tocitem">
    <li id="{@target}">
      <xsl:attribute name="class">
        <xsl:choose>
          <xsl:when test="parent::toc">toc_expanded</xsl:when>
          <xsl:when test="tocitem">toc_collapsed</xsl:when>
        </xsl:choose>
      </xsl:attribute>
      <div class="node_element">
        <div class="wrapper img_wrapper"> </div>
        <div class="wrapper text_wrapper">
          <xsl:variable name="itemurl">
            <xsl:choose>
              <xsl:when test="starts-with(@target,'http://')">
                <xsl:value-of select="@target"/>
              </xsl:when>
              <xsl:when test="starts-with(@target,'https://')">
                <xsl:value-of select="@target"/>
              </xsl:when>
              <xsl:otherwise>
                <xsl:value-of select="concat($helpdir,'/',@target)"/>
              </xsl:otherwise>
            </xsl:choose>
          </xsl:variable>
          <a href="{$itemurl}" class="corrected_url">
            <xsl:value-of select="normalize-space(text())"/>
            <xsl:if test="parent::toc">
              <xsl:text> (Supplemental Software)</xsl:text>
            </xsl:if>
          </a>
        </div>
      </div>
      <xsl:if test="tocitem">
        <ul>
          <xsl:apply-templates select="tocitem"/>
        </ul>
      </xsl:if>
    </li>
  </xsl:template>

</xsl:stylesheet>