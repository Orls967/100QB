# frozen_string_literal: true

require 'json'
require 'fileutils'

module FactorPortfolioLab
  # Deterministic JSON output serializer.
  module JsonOut
    module_function

    def write(path, data)
      FileUtils.mkdir_p(File.dirname(path))
      File.write(path, "#{JSON.pretty_generate(sanitize(data))}\n")
    end

    def sanitize(obj)
      case obj
      when Hash
        obj.each_with_object({}) { |(k, v), h| h[k.to_s] = sanitize(v) }
      when Array
        obj.map { |v| sanitize(v) }
      when Float
        obj.finite? ? obj : nil
      else
        obj
      end
    end
  end
end
